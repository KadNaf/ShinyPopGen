// ============================================================================
// fst_permute_boostrap_OpenMP.cpp
// R-free OpenMP kernels + thin Rcpp wrappers.
// ============================================================================

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <random>
#include <cstring>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ============================================================================
// Structs
// ============================================================================
struct WC84Comp {
  double a, b, c;
  bool ok;
};

struct WC84LocusStats {
  double a, b, c;
  double Ho, Hs, Ht;
  bool ok;
};

// ============================================================================
// Genotype codec (base encoding)
// ============================================================================
inline bool decode_gt_base(int gt, int& a1, int& a2, int base) {
  if (gt == NA_INTEGER || gt <= 0) return false;
  a1 = gt / base;
  a2 = gt % base;
  return (a1 > 0 && a2 > 0);
}

inline bool decode_gt_base_with_missing(int gt, int missing_code, int& a1, int& a2, int base) {
  if (gt == missing_code) return false;
  return decode_gt_base(gt, a1, a2, base);
}

inline int encode_gt_base(int a1, int a2, int base) {
  return a1 * base + a2;
}

// ============================================================================
// Column-major matrix access (R matrices are column-major)
// ============================================================================
static inline int mat_get_colmajor(const int* ptr, int ld, int i, int j) {
  return ptr[i + (size_t)ld * (size_t)j];
}
static inline void mat_set_colmajor(int* ptr, int ld, int i, int j, int v) {
  ptr[i + (size_t)ld * (size_t)j] = v;
}

// ============================================================================
// NA helpers
// ============================================================================
static inline bool is_na_double(double v) {
  return std::isnan(v); // NA_REAL is NaN payload
}

static inline double mean_across_loci_ptr(const double* x, int L) {
  double s = 0.0; int k = 0;
  for (int i = 0; i < L; ++i) {
    const double v = x[i];
    if (is_na_double(v)) continue;
    s += v; k++;
  }
  return (k > 0) ? (s / (double)k) : NA_REAL;
}

// ============================================================================
// Strata builder (pop -> row indices) for "within_pop_alleles" permutations
// ============================================================================
static std::vector<std::vector<int>> build_pop_rows(const IntegerMatrix& mat, int pop_col) {
  const int n = mat.nrow();
  std::unordered_map<int, std::vector<int>> pop_map;
  pop_map.reserve(64);
  
  for (int i = 0; i < n; ++i) pop_map[ mat(i, pop_col) ].push_back(i);
  
  std::vector<std::vector<int>> rows;
  rows.reserve(pop_map.size());
  for (auto& kv : pop_map) rows.push_back(std::move(kv.second));
  return rows;
}

// ============================================================================
// R-FREE permutation helpers operating on raw int* + std::mt19937_64
// ============================================================================
static inline void permute_pop_labels_once_rng_buf(
    int* dat_ptr, int n, int p, int pop_col, std::mt19937_64& rng) {
  
  (void)p;
  std::vector<int> labels((size_t)n);
  for (int i = 0; i < n; ++i) labels[i] = mat_get_colmajor(dat_ptr, n, i, pop_col);
  std::shuffle(labels.begin(), labels.end(), rng);
  for (int i = 0; i < n; ++i) mat_set_colmajor(dat_ptr, n, i, pop_col, labels[i]);
}

static inline void permute_within_pops_once_rng_buf(
    int* dat_ptr, int n, int p, int pop_col,
    int missing_code, int base,
    const std::vector<std::vector<int>>& strata,
    std::mt19937_64& rng) {
  
  std::vector<int> alleles;
  std::vector<int> valid_rows;
  alleles.reserve(256);
  valid_rows.reserve(128);
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    for (const auto& rows : strata) {
      if (rows.size() < 2) continue;
      
      alleles.clear();
      valid_rows.clear();
      
      for (int idx : rows) {
        int a1, a2;
        const int gt = mat_get_colmajor(dat_ptr, n, idx, j);
        if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
        alleles.push_back(a1);
        alleles.push_back(a2);
        valid_rows.push_back(idx);
      }
      
      if (alleles.size() < 2 || valid_rows.empty()) continue;
      std::shuffle(alleles.begin(), alleles.end(), rng);
      
      int k = 0;
      const int aN = (int)alleles.size();
      for (int idx : valid_rows) {
        if (k + 1 >= aN) break;
        mat_set_colmajor(dat_ptr, n, idx, j, encode_gt_base(alleles[k], alleles[k + 1], base));
        k += 2;
      }
    }
  }
}

// ============================================================================
// Pointer-based WC84 locus stats (a,b,c + Ho/Hs/Ht)
// ============================================================================
static inline WC84LocusStats wc84_locus_stats_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int locus_col, int pop_col, int missing_code, int base
) {
  (void)p;
  
  std::unordered_map<int,int> pop_to_idx;
  pop_to_idx.reserve(64);
  
  std::vector<int> pop_labels;
  std::vector<int> n_i;
  std::vector< std::unordered_map<int,int> > acount;
  std::vector< std::unordered_map<int,int> > mho;
  
  auto get_pop_index = [&](int pop)->int {
    auto it = pop_to_idx.find(pop);
    if (it != pop_to_idx.end()) return it->second;
    int idx = (int)pop_labels.size();
    pop_to_idx[pop] = idx;
    pop_labels.push_back(pop);
    n_i.push_back(0);
    acount.emplace_back();
    mho.emplace_back();
    return idx;
  };
  
  std::unordered_map<int,char> all_alleles;
  all_alleles.reserve(128);
  
  int total_n = 0;
  int total_het = 0;
  std::unordered_map<int,int> pooled_acount;
  pooled_acount.reserve(256);
  
  for (int i = 0; i < n_used; ++i) {
    int a1, a2;
    const int gt = mat_get_colmajor(dat_ptr, ld, i, locus_col);
    if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
    
    const int pop = mat_get_colmajor(dat_ptr, ld, i, pop_col);
    const int pi  = get_pop_index(pop);
    
    n_i[pi] += 1;
    total_n += 1;
    
    acount[pi][a1] += 1; acount[pi][a2] += 1;
    pooled_acount[a1] += 1; pooled_acount[a2] += 1;
    all_alleles[a1] = 1; all_alleles[a2] = 1;
    
    if (a1 != a2) {
      total_het += 1;
      mho[pi][a1] += 1;
      mho[pi][a2] += 1;
    }
  }
  
  const int r = (int)pop_labels.size();
  if (r < 2) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  double n_total = 0.0, sum_n2 = 0.0;
  for (int i = 0; i < r; ++i) { n_total += n_i[i]; sum_n2 += (double)n_i[i] * (double)n_i[i]; }
  if (n_total <= 1.0) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  const double n_bar = n_total / (double)r;
  const double n_c   = (n_total - (sum_n2 / n_total)) / (double)(r - 1);
  if (n_c <= 0.0 || n_bar <= 1.0) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  double a_sum = 0.0, b_sum = 0.0, c_sum = 0.0;
  
  for (const auto& kv : all_alleles) {
    const int allele = kv.first;
    
    double pbar_num = 0.0;
    for (int i = 0; i < r; ++i) {
      if (n_i[i] <= 0) continue;
      double pi = 0.0;
      auto it = acount[i].find(allele);
      if (it != acount[i].end()) pi = (double)it->second / (2.0 * (double)n_i[i]);
      pbar_num += (double)n_i[i] * pi;
    }
    const double pbar = pbar_num / n_total;
    
    double s2_num = 0.0;
    for (int i = 0; i < r; ++i) {
      if (n_i[i] <= 0) continue;
      double pi = 0.0;
      auto it = acount[i].find(allele);
      if (it != acount[i].end()) pi = (double)it->second / (2.0 * (double)n_i[i]);
      s2_num += (double)n_i[i] * (pi - pbar) * (pi - pbar);
    }
    const double s2 = s2_num / ((double)(r - 1) * n_bar);
    
    double hbar = 0.0;
    for (int i = 0; i < r; ++i) {
      auto it = mho[i].find(allele);
      if (it != mho[i].end()) hbar += (double)it->second;
    }
    hbar /= n_total;
    
    const double term = pbar * (1.0 - pbar);
    
    const double a = (n_bar / n_c) *
      (s2 - (1.0 / (n_bar - 1.0)) *
      (term - ((double)(r - 1) / (double)r) * s2 - 0.25 * hbar));
    
    const double b = (n_bar / (n_bar - 1.0)) *
      (term - ((double)(r - 1) / (double)r) * s2 - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * hbar);
    
    const double c = 0.5 * hbar;
    
    a_sum += a; b_sum += b; c_sum += c;
  }
  
  const double denom = a_sum + b_sum + c_sum;
  if (!(denom > 0.0)) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  const double Ho = (total_n > 0) ? ((double)total_het / (double)total_n) : NA_REAL;
  
  // Nei HS: unweighted mean of per-population unbiased gene diversities.
  // Each population contributes equally regardless of sample size.
  double Hs_num = 0.0, Hs_den = 0.0;
  for (int pi = 0; pi < r; ++pi) {
    const int ni = n_i[pi];
    if (ni <= 1) continue;  // need n > 1 for unbiased estimate
    const double denom2 = 2.0 * (double)ni;
    double sum_p2 = 0.0;
    for (auto &kv2 : acount[pi]) {
      const double pk = (double)kv2.second / denom2;
      sum_p2 += pk * pk;
    }
    double Hpop = 1.0 - sum_p2;
    Hpop *= (denom2 / (denom2 - 1.0));  // unbiased correction
    Hs_num += Hpop;  // equal weight per population (NOT n_i-weighted)
    Hs_den += 1.0;
  }
  const double Hs = (Hs_den > 0.0) ? (Hs_num / Hs_den) : NA_REAL;
  
  // Nei HT: based on arithmetic mean of per-population allele frequencies.
  // Each population contributes equally regardless of sample size.
  double Ht = NA_REAL;
  {
    int r_active = 0;
    for (int pi = 0; pi < r; ++pi) if (n_i[pi] > 0) ++r_active;
    
    if (r_active > 0 && total_n > 0) {
      // arithmetic mean allele frequency: p̄_a = (1/r) * Σ_i p_{i,a}
      std::unordered_map<int,double> pbar;
      pbar.reserve(128);
      for (int pi = 0; pi < r; ++pi) {
        if (n_i[pi] <= 0) continue;
        const double denom2 = 2.0 * (double)n_i[pi];
        for (auto &kv : acount[pi]) {
          pbar[kv.first] += (double)kv.second / denom2;
        }
      }
      double sum_p2 = 0.0;
      for (auto &kv : pbar) {
        const double p = kv.second / (double)r_active;
        sum_p2 += p * p;
      }
      Ht = 1.0 - sum_p2;
      // unbiased correction using total n across all individuals
      const double denom_tot = 2.0 * (double)total_n;
      if (denom_tot > 1.0) Ht *= (denom_tot / (denom_tot - 1.0));
    }
  }
  return {a_sum, b_sum, c_sum, Ho, Hs, Ht, true};
}

// ============================================================================
// Derived pointer-based helpers (FST per locus + ratio-of-sums overall)
// ============================================================================
static inline void wc84_fst_hs_ht_all_loci_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int pop_col, int missing_code, int base,
    double* fst_out, double* hs_out, double* ht_out
) {
  int out_col = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat_ptr, n_used, ld, p, j, pop_col, missing_code, base);
    if (!s.ok) {
      fst_out[out_col] = NA_REAL;
      hs_out[out_col]  = NA_REAL;
      ht_out[out_col]  = NA_REAL;
    } else {
      const double denom = s.a + s.b + s.c;
      fst_out[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
      hs_out[out_col]  = s.Hs;
      ht_out[out_col]  = s.Ht;
    }
    out_col++;
  }
}

static inline double fst_overall_ratio_of_sums_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int pop_col, int missing_code, int base
) {
  double A = 0.0, D = 0.0;
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat_ptr, n_used, ld, p, j, pop_col, missing_code, base);
    if (!s.ok) continue;
    
    const double denom = s.a + s.b + s.c;
    if (!(denom > 0.0)) continue;
    
    A += s.a;
    D += denom;
  }
  
  return (D > 0.0) ? (A / D) : NA_REAL;
}

// ============================================================================
// G statistic (log-likelihood ratio) for differentiation — FSTAT's test stat
// G = 2 * Σ_pop Σ_allele O_{j,a} * ln(O_{j,a} / E_{j,a})
// where E_{j,a} = 2*n_j * (total_a / 2*N)  (expected under panmixia)
// Returns sum of G across all loci (the global test statistic).
// ============================================================================
static inline double g_stat_differentiation_locus_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int locus_col, int pop_col, int missing_code, int base
) {
  std::unordered_map<int,int> pop_to_idx;
  pop_to_idx.reserve(32);
  std::vector<int> pop_labels;
  std::vector<int> n_i;
  std::vector<std::unordered_map<int,int>> acnt;
  
  auto get_pop_idx = [&](int pop)->int {
    auto it = pop_to_idx.find(pop);
    if (it != pop_to_idx.end()) return it->second;
    int idx = (int)pop_labels.size();
    pop_to_idx[pop] = idx;
    pop_labels.push_back(pop);
    n_i.push_back(0);
    acnt.emplace_back();
    return idx;
  };
  
  std::unordered_map<int,int> total_acnt;
  total_acnt.reserve(64);
  int total_n = 0;
  
  for (int i = 0; i < n_used; ++i) {
    int a1, a2;
    const int gt = mat_get_colmajor(dat_ptr, ld, i, locus_col);
    if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
    const int pop = mat_get_colmajor(dat_ptr, ld, i, pop_col);
    const int pi  = get_pop_idx(pop);
    ++n_i[pi]; ++total_n;
    ++acnt[pi][a1]; ++acnt[pi][a2];
    ++total_acnt[a1]; ++total_acnt[a2];
  }
  
  const int r = (int)pop_labels.size();
  if (r < 2 || total_n < 2) return NA_REAL;
  
  const double two_N = 2.0 * (double)total_n;
  double G = 0.0;
  
  for (int pi = 0; pi < r; ++pi) {
    if (n_i[pi] <= 0) continue;
    const double two_nj = 2.0 * (double)n_i[pi];
    for (auto &kv : acnt[pi]) {
      const int allele = kv.first;
      const double O = (double)kv.second;
      if (O <= 0.0) continue;
      auto it = total_acnt.find(allele);
      if (it == total_acnt.end()) continue;
      const double E = two_nj * ((double)it->second / two_N);
      if (E <= 0.0) continue;
      G += O * std::log(O / E);
    }
  }
  return 2.0 * G;
}

static inline double g_stat_all_loci_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int pop_col, int missing_code, int base
) {
  double G_total = 0.0;
  bool any_ok = false;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    const double g = g_stat_differentiation_locus_ptr_ld(
      dat_ptr, n_used, ld, p, j, pop_col, missing_code, base);
    if (std::isfinite(g)) { G_total += g; any_ok = true; }
  }
  return any_ok ? G_total : NA_REAL;
}


static inline void locus_names_from_colnames(const IntegerMatrix& dat, int pop_col, CharacterVector& out);

// [[Rcpp::export]]
Rcpp::List nei_het_stats_cpp(const Rcpp::IntegerMatrix& dat,
                             int pop_col_1based = 1,
                             int missing_code = 0,
                             int base = 1000) {
  
  if (dat.ncol() < 2)
    Rcpp::stop("Need at least 2 columns: pop + >=1 locus.");
  
  if (base <= 1)
    Rcpp::stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  
  if (pop_col < 0 || pop_col >= dat.ncol())
    Rcpp::stop("pop_col out of range.");
  
  const int N = dat.nrow();
  const int p = dat.ncol();
  const int L = p - 1;
  
  Rcpp::NumericVector Ho(L, NA_REAL);
  Rcpp::NumericVector Hs(L, NA_REAL);
  Rcpp::NumericVector Ht(L, NA_REAL);
  
  Rcpp::IntegerVector k_eff(L);
  std::fill(k_eff.begin(), k_eff.end(), NA_INTEGER);
  
  // map pop label -> index
  std::unordered_map<int,int> pop_to_idx;
  pop_to_idx.reserve(64);
  
  std::vector<int> pop_labels;
  pop_labels.reserve(16);
  
  auto get_pop_index = [&](int pop)->int {
    auto it = pop_to_idx.find(pop);
    if (it != pop_to_idx.end()) return it->second;
    
    int idx = (int)pop_labels.size();
    pop_to_idx[pop] = idx;
    pop_labels.push_back(pop);
    return idx;
  };
  
  // pre-scan populations
  for (int i = 0; i < N; ++i) {
    int pop = dat(i, pop_col);
    get_pop_index(pop);
  }
  
  const int r = (int)pop_labels.size();
  
  if (r < 1)
    return Rcpp::List::create(_["Ho"]=Ho, _["Hs"]=Hs, _["Ht"]=Ht);
  
  int out_col = 0;
  
  for (int j = 0; j < p; ++j) {
    
    if (j == pop_col) continue;
    
    std::vector< std::unordered_map<int,int> > acount(r);
    std::vector<int> n_i(r, 0);
    std::vector<int> het_i(r, 0);
    
    std::unordered_map<int,int> total_acount;
    total_acount.reserve(256);
    
    int total_n = 0;
    int total_het = 0;
    
    for (int i = 0; i < N; ++i) {
      
      int a1, a2;
      int gt = dat(i, j);
      
      if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base))
        continue;
      
      int pop = dat(i, pop_col);
      int pi = pop_to_idx[pop];
      
      n_i[pi] += 1;
      total_n += 1;
      
      acount[pi][a1] += 1;
      acount[pi][a2] += 1;
      
      total_acount[a1] += 1;
      total_acount[a2] += 1;
      
      if (a1 != a2) {
        het_i[pi] += 1;
        total_het += 1;
      }
    }
    
    if (total_n == 0) {
      out_col++;
      continue;
    }
    
    int k_locus = 0;
    for (int pi = 0; pi < r; ++pi) {
      if (n_i[pi] > 0) k_locus++;
    }
    
    k_eff[out_col] = k_locus;
    
    // Ho = mean heterozygosity across populations
    double Ho_sum = 0.0;
    double Ho_den = 0.0;
    
    for (int pi = 0; pi < r; ++pi) {
      if (n_i[pi] <= 0) continue;
      
      Ho_sum += (double)het_i[pi] / (double)n_i[pi];
      Ho_den += 1.0;
    }
    
    Ho[out_col] = (Ho_den > 0.0) ? (Ho_sum / Ho_den) : NA_REAL;
    
    // Hs
    double Hs_num = 0.0;
    double Hs_den = 0.0;
    
    for (int pi = 0; pi < r; ++pi) {
      
      if (n_i[pi] <= 1) continue;
      
      const double denom = 2.0 * (double)n_i[pi];
      
      double sum_p2 = 0.0;
      
      for (auto &kv : acount[pi]) {
        
        const double pk = (double)kv.second / denom;
        sum_p2 += pk * pk;
      }
      
      double Hpop = 1.0 - sum_p2;
      
      // unbiased correction
      Hpop *= (denom / (denom - 1.0));
      
      Hs_num += Hpop;
      Hs_den += 1.0;
    }
    
    if (Hs_den > 0.0)
      Hs[out_col] = Hs_num / Hs_den;
    
    // Ht
      {
        int k_eff = 0;
        
        for (int pi = 0; pi < r; ++pi)
          if (n_i[pi] > 0) k_eff++;
          
          if (k_eff <= 0) {
            
            Ht[out_col] = NA_REAL;
            
          } else {
            
            std::unordered_map<int,double> p_sum;
            p_sum.reserve(256);
            
            for (int pi = 0; pi < r; ++pi) {
              
              if (n_i[pi] <= 0) continue;
              
              const double denom_pi = 2.0 * (double)n_i[pi];
              
              for (auto &kv : acount[pi]) {
                
                const int allele = kv.first;
                const double p_aj = (double)kv.second / denom_pi;
                
                p_sum[allele] += p_aj;
              }
            }
            
            double sum_p2 = 0.0;
            
            for (auto &kv : p_sum) {
              
              const double pbar = kv.second / (double)k_eff;
              sum_p2 += pbar * pbar;
            }
            
            double Htot = 1.0 - sum_p2;
            
            const double denom_tot = 2.0 * (double)total_n;
            
            if (denom_tot > 1.0)
              Htot *= (denom_tot / (denom_tot - 1.0));
            
            Ht[out_col] = Htot;
          }
      }
      
      out_col++;
  }
  
  // attach locus names
  Rcpp::CharacterVector cn = colnames(dat);
  Rcpp::CharacterVector loc_names(L);
  
  int kk = 0;
  const bool has_colnames = (cn.size() == p);
  
  for (int j = 0; j < p; ++j) {
    
    if (j == pop_col) continue;
    
    if (has_colnames)
      loc_names[kk] = cn[j];
    else
      loc_names[kk] = std::string("L") + std::to_string(kk + 1);
    
    kk++;
  }
  
  Ho.attr("names") = loc_names;
  Hs.attr("names") = loc_names;
  Ht.attr("names") = loc_names;
  
  return Rcpp::List::create(
    _["Ho"] = Ho,
    _["Hs"] = Hs,
    _["Ht"] = Ht,
    _["k_eff"] = k_eff,
    _["locus_names"] = loc_names
  );
}

// ============================================================================
// R-FREE OpenMP permutation kernel (FST only)
// ============================================================================
static void fst_perm_kernel_parallel_rfree(
    const int* dat_ptr, int n, int p,
    int pop_col, int missing_code, int base,
    int perm_mode,
    int B, int n_threads, uint64_t seed0,
    const std::vector<std::vector<int>>& strata,
    std::vector<double>& fst_perm_buf,
    std::vector<double>& fst_overall_buf,
    std::vector<double>& g_overall_buf,  // B — global G per permutation
    std::vector<double>& g_perm_buf      // B*L — per-locus G per permutation
) {
  const int L = p - 1;

  int T = std::max(1, n_threads);
#ifdef _OPENMP
  omp_set_num_threads(T);
#else
  T = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<int> pm((size_t)n * (size_t)p);
  std::vector<double> fst_tmp((size_t)L, NA_REAL);
  std::vector<double> g_tmp((size_t)L, NA_REAL);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    std::mt19937_64 rng(sb);

    std::memcpy(pm.data(), dat_ptr, sizeof(int) * (size_t)n * (size_t)p);

    if (perm_mode == 1) {
      permute_pop_labels_once_rng_buf(pm.data(), n, p, pop_col, rng);
    } else {
      permute_within_pops_once_rng_buf(pm.data(), n, p, pop_col, missing_code, base, strata, rng);
    }

    int out_col = 0;
    for (int j = 0; j < p; ++j) {
      if (j == pop_col) continue;
      WC84LocusStats s = wc84_locus_stats_ptr_ld(pm.data(), n, n, p, j, pop_col, missing_code, base);
      if (!s.ok) {
        fst_tmp[out_col] = NA_REAL;
        g_tmp[out_col]   = NA_REAL;
      } else {
        const double denom = s.a + s.b + s.c;
        fst_tmp[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
        g_tmp[out_col]   = g_stat_differentiation_locus_ptr_ld(pm.data(), n, n, p, j, pop_col, missing_code, base);
      }
      out_col++;
    }

    const double fst_overall = fst_overall_ratio_of_sums_ptr_ld(pm.data(), n, n, p, pop_col, missing_code, base);
    const double g_overall = g_stat_all_loci_ptr_ld(pm.data(), n, n, p, pop_col, missing_code, base);
    g_overall_buf[(size_t)b] = g_overall;
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) {
      fst_perm_buf[row0 + (size_t)ell] = fst_tmp[(size_t)ell];
      g_perm_buf[row0 + (size_t)ell]   = g_tmp[(size_t)ell];
    }
    fst_overall_buf[(size_t)b] = fst_overall;
  }
}
}

// ============================================================================
// R-FREE OpenMP pop-block bootstrap kernel (FST + HS + HT)
// ============================================================================
static void boot_popblock_kernel_parallel_rfree(
    const int* dat_ptr, int n, int p, int pop_col,
    int missing_code, int base,
    int B, int n_threads, uint64_t seed0,
    const std::vector<std::vector<int>>& base_rows, // population blocks (indices in original)
    std::vector<double>& fst_boot_buf,              // B*L
    std::vector<double>& fst_overall_buf,           // B
    std::vector<double>& hs_boot_buf,               // B*L
    std::vector<double>& hs_overall_buf,            // B
    std::vector<double>& ht_boot_buf,               // B*L
    std::vector<double>& ht_overall_buf             // B
) {
  const int L = p - 1;
  const int P = (int)base_rows.size();
  if (P <= 0) return;
  
  int max_block = 0;
  for (const auto& rows : base_rows) max_block = std::max(max_block, (int)rows.size());
  const int ld = std::max(1, max_block * P);
  
  int T = std::max(1, n_threads);
#ifdef _OPENMP
  omp_set_num_threads(T);
#else
  T = 1;
#endif
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<int> mbuf((size_t)ld * (size_t)p);
  std::vector<double> fst_tmp((size_t)L), hs_tmp((size_t)L), ht_tmp((size_t)L);
  
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    std::mt19937_64 rng(sb);
    std::uniform_int_distribution<int> U(0, P - 1);
    
    int rr = 0;
    for (int k = 0; k < P; ++k) {
      const int pick = U(rng);
      const auto& rows = base_rows[(size_t)pick];
      const int block_id = k + 1;
      
      for (int rid : rows) {
        for (int j = 0; j < p; ++j) {
          const int v = mat_get_colmajor(dat_ptr, n, rid, j);
          mat_set_colmajor(mbuf.data(), ld, rr, j, v);
        }
        mat_set_colmajor(mbuf.data(), ld, rr, pop_col, block_id);
        rr++;
        if (rr >= ld) break;
      }
      if (rr >= ld) break;
    }
    const int n_used = rr;
    
    wc84_fst_hs_ht_all_loci_ptr_ld(
      mbuf.data(), n_used, ld, p,
      pop_col, missing_code, base,
      fst_tmp.data(), hs_tmp.data(), ht_tmp.data()
    );
    
    const double fst_overall = fst_overall_ratio_of_sums_ptr_ld(mbuf.data(), n_used, ld, p, pop_col, missing_code, base);
    const double hs_overall  = mean_across_loci_ptr(hs_tmp.data(), L);
    const double ht_overall  = mean_across_loci_ptr(ht_tmp.data(), L);
    
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) {
      fst_boot_buf[row0 + (size_t)ell] = fst_tmp[(size_t)ell];
      hs_boot_buf [row0 + (size_t)ell] = hs_tmp [(size_t)ell];
      ht_boot_buf [row0 + (size_t)ell] = ht_tmp [(size_t)ell];
    }
    fst_overall_buf[(size_t)b] = fst_overall;
    hs_overall_buf [(size_t)b] = hs_overall;
    ht_overall_buf [(size_t)b] = ht_overall;
  }
}
}

// ============================================================================
// Export: observed WC84 stats (per-locus + overall) using pointer core
// ============================================================================
static inline void locus_names_from_colnames(const IntegerMatrix& dat, int pop_col, CharacterVector& out) {
  const int p = dat.ncol();
  const int L = p - 1;
  CharacterVector cn = colnames(dat);
  
  out = CharacterVector(L);
  int k = 0;
  const bool has = (cn.size() == p);
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    if (has) {
      out[k] = cn[j];                        // same type: Rcpp::String proxy
    } else {
      out[k] = Rcpp::String(std::string("L") + std::to_string(k + 1));
    }
    ++k;
  }
}

// [[Rcpp::export]]
List observed_wc84_stats_cpp(const IntegerMatrix& dat,
                             int pop_col_1based = 1,
                             int missing_code = 0,
                             int base = 1000) {
  
  if (dat.ncol() < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (base <= 1) stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= dat.ncol()) stop("pop_col out of range.");
  
  const int n = dat.nrow();
  const int p = dat.ncol();
  const int L = p - 1;
  
  NumericVector FST(L, NA_REAL), FIT(L, NA_REAL), FIS(L, NA_REAL);
  NumericVector Ho(L, NA_REAL), HS(L, NA_REAL), HT(L, NA_REAL);
  
  double A_sum = 0.0, B_sum = 0.0, C_sum = 0.0;
  
  int out_col = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat.begin(), n, n, p, j, pop_col, missing_code, base);
    if (s.ok) {
      const double denom = s.a + s.b + s.c;
      const double denom_fis = s.b + s.c;
      
      FST[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
      FIT[out_col] = (denom > 0.0) ? ((s.a + s.b) / denom) : NA_REAL;
      FIS[out_col] = (denom_fis > 0.0) ? (s.b / denom_fis) : NA_REAL;
      
      Ho[out_col]  = s.Ho;
      HS[out_col]  = s.Hs;
      HT[out_col]  = s.Ht;
      
      if (denom > 0.0) {
        A_sum += s.a;
        B_sum += s.b;
        C_sum += s.c;
      }
    }
    out_col++;
  }
  
  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);
  
  FST.attr("names") = loc_names;
  FIT.attr("names") = loc_names;
  FIS.attr("names") = loc_names;
  Ho.attr("names")  = loc_names;
  HS.attr("names")  = loc_names;
  HT.attr("names")  = loc_names;
  
  const double denom_overall = A_sum + B_sum + C_sum;
  const double denom_fis_overall = B_sum + C_sum;
  
  const double FST_overall = (denom_overall > 0.0) ? (A_sum / denom_overall) : NA_REAL;
  const double FIT_overall = (denom_overall > 0.0) ? ((A_sum + B_sum) / denom_overall) : NA_REAL;
  const double FIS_overall = (denom_fis_overall > 0.0) ? (B_sum / denom_fis_overall) : NA_REAL;
  
  return List::create(
    _["FST"] = FST,
    _["FIT"] = FIT,
    _["FIS"] = FIS,
    _["HI"]  = Ho,
    _["HS"]  = HS,
    _["HT"]  = HT,
    _["FST_overall_ratio_of_sums"] = FST_overall,
    _["FIT_overall_ratio_of_sums"] = FIT_overall,
    _["FIS_overall_ratio_of_sums"] = FIS_overall,
    _["locus_names"] = loc_names
  );
}

// ============================================================================
// Export: per-locus WC84 variance components for locus bootstrap
// Returns a data frame with columns A, B, C, HS, HT per locus.
// The caller resamples rows (loci) with replacement and recomputes
// ratio-of-sums to obtain the locus-bootstrap distribution.
// ============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame wc84_locus_components_cpp(
    const Rcpp::IntegerMatrix& dat,
    int pop_col_1based = 1,
    int missing_code = 0,
    int base = 1000
) {
  if (dat.ncol() < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (base <= 1) stop("base must be > 1.");

  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= dat.ncol()) stop("pop_col out of range.");

  const int n = dat.nrow();
  const int p = dat.ncol();
  const int L = p - 1;

  NumericVector A_v(L, NA_REAL), B_v(L, NA_REAL), C_v(L, NA_REAL);
  NumericVector HS_v(L, NA_REAL), HT_v(L, NA_REAL);

  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);

  int out = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    WC84LocusStats s = wc84_locus_stats_ptr_ld(
      dat.begin(), n, n, p, j, pop_col, missing_code, base);
    if (s.ok) {
      A_v[out]  = s.a;
      B_v[out]  = s.b;
      C_v[out]  = s.c;
      HS_v[out] = s.Hs;
      HT_v[out] = s.Ht;
    }
    ++out;
  }

  return DataFrame::create(
    _["locus"] = loc_names,
    _["A"]     = A_v,
    _["B"]     = B_v,
    _["C"]     = C_v,
    _["HS"]    = HS_v,
    _["HT"]    = HT_v
  );
}

// ============================================================================
// Export: parallel permutation (FST)
// ============================================================================
static inline uint64_t seed0_from_double(double seed) {
  const uint64_t s = (uint64_t)(seed < 0 ? -seed : seed);
  return s + 0x9e3779b97f4a7c15ULL;
}

// ============================================================================
// Locus bootstrap for global WC84 estimators (OpenMP).
// Input: per-locus WC84 variance components (A, B, C) and gene diversities
//        (HS, HT), already computed by wc84_locus_components_cpp().
// Each replicate resamples L loci with replacement and recomputes the
// ratio-of-sums: FST = ΣA/Σ(A+B+C), FIT = Σ(A+B)/Σ(A+B+C),
// FIS = ΣB/Σ(B+C), HS = mean(HS_l), HT = mean(HT_l).
// Returns a 5-row summary DataFrame (one row per statistic).
// ============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame locus_bootstrap_wc84_cpp(
    const Rcpp::NumericVector& A,
    const Rcpp::NumericVector& Bv,
    const Rcpp::NumericVector& C,
    const Rcpp::NumericVector& HS,
    const Rcpp::NumericVector& HT,
    int    B_reps     = 1000,
    double conf_level = 0.95,
    double seed       = 1.0,
    int    n_threads  = 1
) {
  const int L_all = A.size();
  if (L_all == 0) stop("Empty component vectors.");
  if (Bv.size() != L_all || C.size() != L_all ||
      HS.size() != L_all || HT.size() != L_all)
    stop("All component vectors must have the same length.");
  if (B_reps <= 0)                          stop("B_reps must be positive.");
  if (conf_level <= 0.0 || conf_level >= 1.0)
    stop("conf_level must be in (0, 1).");

  // Filter to loci with all five components finite
  std::vector<double> a_v, b_v, c_v, hs_v, ht_v;
  a_v.reserve((size_t)L_all); b_v.reserve((size_t)L_all);
  c_v.reserve((size_t)L_all); hs_v.reserve((size_t)L_all);
  ht_v.reserve((size_t)L_all);
  for (int i = 0; i < L_all; ++i) {
    if (std::isfinite(A[i])  && std::isfinite(Bv[i]) &&
        std::isfinite(C[i])  && std::isfinite(HS[i]) &&
        std::isfinite(HT[i])) {
      a_v.push_back(A[i]);  b_v.push_back(Bv[i]); c_v.push_back(C[i]);
      hs_v.push_back(HS[i]); ht_v.push_back(HT[i]);
    }
  }
  const int L = (int)a_v.size();
  if (L < 2) stop("Need at least 2 valid loci for locus bootstrap.");

  // Observed ratio-of-sums from the full valid locus set
  double sum_a = 0.0, sum_b = 0.0, sum_c = 0.0,
         sum_hs = 0.0, sum_ht = 0.0;
  for (int i = 0; i < L; ++i) {
    sum_a  += a_v[(size_t)i];  sum_b  += b_v[(size_t)i];
    sum_c  += c_v[(size_t)i];  sum_hs += hs_v[(size_t)i];
    sum_ht += ht_v[(size_t)i];
  }
  const double dABC   = sum_a + sum_b + sum_c;
  const double dBC    = sum_b + sum_c;
  const double Ld     = (double)L;
  const double obs_fst = (dABC > 0.0) ? sum_a / dABC        : NA_REAL;
  const double obs_fit = (dABC > 0.0) ? (sum_a+sum_b)/dABC  : NA_REAL;
  const double obs_fis = (dBC  > 0.0) ? sum_b / dBC         : NA_REAL;
  const double obs_hs  = sum_hs / Ld;
  const double obs_ht  = sum_ht / Ld;

  // Bootstrap storage (one value per replicate, per statistic)
  std::vector<double> boot_fst((size_t)B_reps, NA_REAL);
  std::vector<double> boot_fit((size_t)B_reps, NA_REAL);
  std::vector<double> boot_fis((size_t)B_reps, NA_REAL);
  std::vector<double> boot_hs ((size_t)B_reps, NA_REAL);
  std::vector<double> boot_ht ((size_t)B_reps, NA_REAL);

  const uint64_t seed0 = seed0_from_double(seed);

  int T = std::max(1, n_threads);
#ifdef _OPENMP
  omp_set_num_threads(T);
#else
  T = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int b = 0; b < B_reps; ++b) {
      // Per-replicate RNG — same seeding pattern as the rest of this file
      const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
      std::mt19937_64 rng(sb);
      std::uniform_int_distribution<int> U(0, L - 1);

      double sa = 0.0, sb2 = 0.0, sc = 0.0, shs = 0.0, sht = 0.0;
      for (int k = 0; k < L; ++k) {
        const int idx = U(rng);
        sa  += a_v[(size_t)idx];
        sb2 += b_v[(size_t)idx];
        sc  += c_v[(size_t)idx];
        shs += hs_v[(size_t)idx];
        sht += ht_v[(size_t)idx];
      }
      const double dabc = sa + sb2 + sc;
      const double dbc  = sb2 + sc;
      boot_fst[(size_t)b] = (dabc > 0.0) ? sa / dabc          : NA_REAL;
      boot_fit[(size_t)b] = (dabc > 0.0) ? (sa+sb2) / dabc    : NA_REAL;
      boot_fis[(size_t)b] = (dbc  > 0.0) ? sb2 / dbc          : NA_REAL;
      boot_hs [(size_t)b] = shs / Ld;
      boot_ht [(size_t)b] = sht / Ld;
    }
  }

  // Summarise bootstrap distribution -> Observed, Boot_Mean, SE, CI_L, CI_U
  const double alpha = (1.0 - conf_level) / 2.0;

  auto summarise = [&](const std::vector<double>& x, double obs)
      -> std::array<double, 5> {
    std::vector<double> v;
    v.reserve(x.size());
    for (double xi : x) if (std::isfinite(xi)) v.push_back(xi);
    if ((int)v.size() < 2) return {obs, NA_REAL, NA_REAL, NA_REAL, NA_REAL};
    double s = 0.0, s2 = 0.0;
    for (double xi : v) { s += xi; s2 += xi * xi; }
    const double mn = s / (double)v.size();
    const double var = s2/(double)v.size() - mn*mn;
    const double se = std::sqrt(var > 0.0
      ? var * (double)v.size() / (double)(v.size() - 1) : 0.0);
    std::sort(v.begin(), v.end());
    const int nv = (int)v.size();
    auto q = [&](double p) -> double {
      const double fi = p * (double)(nv - 1);
      const int i0 = (int)std::floor(fi);
      const int i1 = std::min(i0 + 1, nv - 1);
      return v[(size_t)i0] * (1.0 - (fi - (double)i0))
           + v[(size_t)i1] * (fi - (double)i0);
    };
    return {obs, mn, se, q(alpha), q(1.0 - alpha)};
  };

  auto r_fst = summarise(boot_fst, obs_fst);
  auto r_fit = summarise(boot_fit, obs_fit);
  auto r_fis = summarise(boot_fis, obs_fis);
  auto r_hs  = summarise(boot_hs,  obs_hs);
  auto r_ht  = summarise(boot_ht,  obs_ht);

  CharacterVector stat_names = {"FST", "FIT", "FIS", "HS", "HT"};
  NumericVector obs_v(5), mean_v(5), se_v(5), ci_l_v(5), ci_u_v(5);
  const std::vector<std::array<double,5>> rows =
    {r_fst, r_fit, r_fis, r_hs, r_ht};
  for (int i = 0; i < 5; ++i) {
    obs_v[i]  = rows[(size_t)i][0];
    mean_v[i] = rows[(size_t)i][1];
    se_v[i]   = rows[(size_t)i][2];
    ci_l_v[i] = rows[(size_t)i][3];
    ci_u_v[i] = rows[(size_t)i][4];
  }

  return DataFrame::create(
    _["Statistic"] = stat_names,
    _["Observed"]  = obs_v,
    _["Boot_Mean"] = mean_v,
    _["SE"]        = se_v,
    _["CI_L"]      = ci_l_v,
    _["CI_U"]      = ci_u_v
  );
}

// [[Rcpp::export]]
List batch_permute_wc84_fst_parallel(const IntegerMatrix& dat,
                                           int pop_col_1based,
                                           int missing_code,
                                           int base,
                                           int B,
                                           int n_threads = 1,
                                           double seed = 1.0,
                                           std::string pval_method = "two_sided_abs",
                                           std::string perm_scheme = "within_pop_alleles") {
  
  const int n = dat.nrow();
  const int p = dat.ncol();
  if (p < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (B <= 0) stop("B must be positive.");
  if (base <= 1) stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= p) stop("pop_col out of range.");
  
  const int L = p - 1;
  
  int perm_mode = 0;
  std::vector<std::vector<int>> strata;
  if (perm_scheme == "within_pop_alleles") {
    perm_mode = 0;
    strata = build_pop_rows(dat, pop_col);
  } else if (perm_scheme == "permute_pop_labels") {
    perm_mode = 1;
  } else {
    stop("Unknown perm_scheme.");
  }
  
  // observed (pointer core)
  NumericVector fst_obs(L, NA_REAL);
  NumericVector g_obs(L, NA_REAL);
  int out_col = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat.begin(), n, n, p, j, pop_col, missing_code, base);
    if (s.ok) {
      const double denom = s.a + s.b + s.c;
      fst_obs[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
      g_obs[out_col]   = g_stat_differentiation_locus_ptr_ld(dat.begin(), n, n, p, j, pop_col, missing_code, base);
    }
    out_col++;
  }
  const double fst_overall_obs = fst_overall_ratio_of_sums_ptr_ld(dat.begin(), n, n, p, pop_col, missing_code, base);
  const double g_overall_obs = g_stat_all_loci_ptr_ld(dat.begin(), n, n, p, pop_col, missing_code, base);
  int T = std::max(1, n_threads);
#ifndef _OPENMP
  T = 1;
#endif
  
  const uint64_t seed0 = seed0_from_double(seed);
  
  std::vector<double> fst_perm_buf((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> fst_overall_buf((size_t)B, NA_REAL);
  std::vector<double> g_overall_buf((size_t)B, NA_REAL);
  std::vector<double> g_perm_buf((size_t)B * (size_t)L, NA_REAL);

  fst_perm_kernel_parallel_rfree(
    dat.begin(), n, p,
    pop_col, missing_code, base,
    perm_mode,
    B, T, seed0,
    strata,
    fst_perm_buf,
    fst_overall_buf,
    g_overall_buf,
    g_perm_buf
  );
  
  // materialize R outputs
  NumericMatrix fst_perm(B, L);
  NumericVector fst_overall_perm(B, NA_REAL);
  for (int b = 0; b < B; ++b) {
    fst_overall_perm[b] = fst_overall_buf[(size_t)b];
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) fst_perm(b, ell) = fst_perm_buf[row0 + (size_t)ell];
  }
  
  // FST-based per-locus p-values (kept for reference)
  NumericVector p_fst(L, NA_REAL);
  for (int ell = 0; ell < L; ++ell) {
    const double obs = fst_obs[ell];
    if (NumericVector::is_na(obs)) continue;

    int count = 0, valid = 0;
    for (int b = 0; b < B; ++b) {
      const double v = fst_perm(b, ell);
      if (NumericVector::is_na(v)) continue;
      valid++;

      if (pval_method == "greater") {
        if (v >= obs) count++;
      } else if (pval_method == "less") {
        if (v <= obs) count++;
      } else {
        if (std::fabs(v) >= std::fabs(obs)) count++;
      }
    }
    if (valid > 0) p_fst[ell] = (1.0 + count) / ((double)valid + 1.0);
  }

  // G-based per-locus p-values (one-sided: G >= 0 by definition)
  NumericMatrix g_perm(B, L);
  for (int b = 0; b < B; ++b) {
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) g_perm(b, ell) = g_perm_buf[row0 + (size_t)ell];
  }

  NumericVector p_G(L, NA_REAL);
  for (int ell = 0; ell < L; ++ell) {
    const double obs_g = g_obs[ell];
    if (NumericVector::is_na(obs_g) || !std::isfinite(obs_g)) continue;
    int count = 0, valid = 0;
    for (int b = 0; b < B; ++b) {
      const double v = g_perm(b, ell);
      if (NumericVector::is_na(v) || !std::isfinite(v)) continue;
      ++valid;
      if (v >= obs_g) ++count;
    }
    if (valid > 0) p_G[ell] = (1.0 + count) / ((double)valid + 1.0);
  }

  // Overall p-value uses G statistic (FSTAT convention), not theta.
  double p_fst_overall = NA_REAL;
  if (std::isfinite(g_overall_obs)) {
    int count = 0, valid = 0;
    for (int b = 0; b < B; ++b) {
      const double v = g_overall_buf[(size_t)b];
      if (!std::isfinite(v)) continue;
      ++valid;
      if (v >= g_overall_obs) ++count;  // one-sided: G_perm >= G_obs
    }
    if (valid > 0) p_fst_overall = (1.0 + count) / ((double)valid + 1.0);
  }
  
  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);
  colnames(fst_perm) = loc_names;
  colnames(g_perm)   = loc_names;
  fst_obs.attr("names") = loc_names;
  g_obs.attr("names")   = loc_names;
  p_fst.attr("names")   = loc_names;
  p_G.attr("names")     = loc_names;
  
  NumericVector g_overall_perm(B, NA_REAL);  //check place
  for (int b = 0; b < B; ++b) g_overall_perm[b] = g_overall_buf[(size_t)b]; //check place
  
  return List::create(
    _["FST_obs"]         = fst_obs,
    _["FST_overall_obs"] = fst_overall_obs,
    _["FST_perm"]        = fst_perm,
    _["FST_overall_perm"]= fst_overall_perm,
    _["G_locus_obs"]     = g_obs,
    _["G_locus_perm"]    = g_perm,
    _["G_overall_obs"]   = g_overall_obs,
    _["G_overall_perm"]  = g_overall_perm,
    _["p_FST"]           = p_fst,
    _["p_G"]             = p_G,
    _["p_FST_overall"]   = p_fst_overall,
    _["locus_names"]     = loc_names,
    _["perm_scheme"]     = perm_scheme,
    _["pval_method"]     = pval_method,
    _["n_perm"]          = B,
    _["n_threads"]       = T,
    _["seed"]            = seed
  );
}

// ============================================================================
// Export: parallel pop-block bootstrap (FST + HS + HT)
// ============================================================================

// [[Rcpp::export]]
List boot_popblock_wc84_parallel(const IntegerMatrix& mat,
                                       int pop_col_1based = 1,
                                       int missing_code = 0,
                                       int base = 1000,
                                       int B = 1000,
                                       int n_threads = 1,
                                       double seed = 1.0) {
  
  const int p = mat.ncol();
  if (p < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (B <= 0) stop("B must be positive.");
  if (base <= 1) stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= p) stop("pop_col out of range.");
  
  const int n = mat.nrow();
  const int L = p - 1;
  
  auto base_rows = build_pop_rows(mat, pop_col);
  
  // observed from pointer core
  std::vector<double> fst_obs_buf((size_t)L, NA_REAL), hs_obs_buf((size_t)L, NA_REAL), ht_obs_buf((size_t)L, NA_REAL);
  wc84_fst_hs_ht_all_loci_ptr_ld(mat.begin(), n, n, p, pop_col, missing_code, base,
                                 fst_obs_buf.data(), hs_obs_buf.data(), ht_obs_buf.data());
  
  NumericVector fst_obs(L), hs_obs(L), ht_obs(L);
  for (int i = 0; i < L; ++i) { fst_obs[i] = fst_obs_buf[i]; hs_obs[i] = hs_obs_buf[i]; ht_obs[i] = ht_obs_buf[i]; }
  
  const double fst_overall_obs = fst_overall_ratio_of_sums_ptr_ld(mat.begin(), n, n, p, pop_col, missing_code, base);
  const double hs_overall_obs  = mean_across_loci_ptr(hs_obs_buf.data(), L);
  const double ht_overall_obs  = mean_across_loci_ptr(ht_obs_buf.data(), L);
  
  int T = std::max(1, n_threads);
#ifndef _OPENMP
  T = 1;
#endif
  const uint64_t seed0 = seed0_from_double(seed);
  
  // buffers
  std::vector<double> fst_boot_buf((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> hs_boot_buf ((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> ht_boot_buf ((size_t)B * (size_t)L, NA_REAL);
  
  std::vector<double> fst_overall_buf((size_t)B, NA_REAL);
  std::vector<double> hs_overall_buf ((size_t)B, NA_REAL);
  std::vector<double> ht_overall_buf ((size_t)B, NA_REAL);
  
  boot_popblock_kernel_parallel_rfree(
    mat.begin(), n, p, pop_col,
    missing_code, base,
    B, T, seed0,
    base_rows,
    fst_boot_buf, fst_overall_buf,
    hs_boot_buf,  hs_overall_buf,
    ht_boot_buf,  ht_overall_buf
  );
  
  // materialize
  NumericMatrix fst_boot(B, L), hs_boot(B, L), ht_boot(B, L);
  NumericVector fst_overall_boot(B), hs_overall_boot(B), ht_overall_boot(B);
  
  for (int b = 0; b < B; ++b) {
    fst_overall_boot[b] = fst_overall_buf[(size_t)b];
    hs_overall_boot[b]  = hs_overall_buf[(size_t)b];
    ht_overall_boot[b]  = ht_overall_buf[(size_t)b];
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) {
      fst_boot(b, ell) = fst_boot_buf[row0 + (size_t)ell];
      hs_boot (b, ell) = hs_boot_buf [row0 + (size_t)ell];
      ht_boot (b, ell) = ht_boot_buf [row0 + (size_t)ell];
    }
  }
  
  CharacterVector loc_names;
  locus_names_from_colnames(mat, pop_col, loc_names);
  colnames(fst_boot) = loc_names;
  colnames(hs_boot)  = loc_names;
  colnames(ht_boot)  = loc_names;
  
  fst_obs.attr("names") = loc_names;
  hs_obs.attr("names")  = loc_names;
  ht_obs.attr("names")  = loc_names;
  
  return List::create(
    _["FST_obs"] = fst_obs,
    _["FST_overall_obs"] = fst_overall_obs,
    _["FST_boot"] = fst_boot,
    _["FST_overall_boot"] = fst_overall_boot,
    
    _["HS_obs"] = hs_obs,
    _["HS_overall_obs"] = hs_overall_obs,
    _["HS_boot"] = hs_boot,
    _["HS_overall_boot"] = hs_overall_boot,
    
    _["HT_obs"] = ht_obs,
    _["HT_overall_obs"] = ht_overall_obs,
    _["HT_boot"] = ht_boot,
    _["HT_overall_boot"] = ht_overall_boot,
    
    _["boot_type"] = "pop_block",
    _["n_boot"] = B,
    _["n_threads"] = T,
    _["seed"] = seed,
    _["locus_names"] = loc_names
  );
}

// ============================================================================
// Individual bootstrap for HS: resample individuals within each population
// Returns HS_boot (B x L) and HS_overall_boot (B) — loci with replacement
// would give locus-level CI, handled separately via locus_bootstrap_wc84_cpp.
// ============================================================================
// [[Rcpp::export]]
List boot_indiv_hs_cpp(
    const IntegerMatrix& dat,
    int pop_col_1based = 1,
    int missing_code   = 0,
    int base           = 1000,
    int B              = 1000,
    double seed        = 1.0,
    int n_threads      = 1
) {
  const int n = dat.nrow();
  const int p = dat.ncol();
  if (p < 2) stop("Need at least 2 columns.");
  if (B <= 0) stop("B must be positive.");
  if (base <= 1) stop("base must be > 1.");

  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= p) stop("pop_col out of range.");
  const int L = p - 1;

  std::vector<std::vector<int>> pop_rows = build_pop_rows(dat, pop_col);
  const int P = (int)pop_rows.size();
  if (P < 1) stop("No populations found.");

  const uint64_t seed0 = seed0_from_double(seed);

  std::vector<double> hs_boot_buf((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> hs_overall_buf((size_t)B, NA_REAL);

  int T = std::max(1, n_threads);
#ifndef _OPENMP
  T = 1;
#endif
#ifdef _OPENMP
  omp_set_num_threads(T);
#pragma omp parallel
#endif
  {
    std::vector<int> mbuf((size_t)n * (size_t)p);
    std::vector<double> fst_tmp((size_t)L, NA_REAL);
    std::vector<double> hs_tmp((size_t)L, NA_REAL);
    std::vector<double> ht_tmp((size_t)L, NA_REAL);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int b = 0; b < B; ++b) {
      const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
      std::mt19937_64 rng(sb);

      int rr = 0;
      for (int k = 0; k < P; ++k) {
        const auto& rows = pop_rows[(size_t)k];
        const int nk = (int)rows.size();
        if (nk == 0) continue;
        std::uniform_int_distribution<int> U_indiv(0, nk - 1);
        for (int i = 0; i < nk; ++i) {
          const int rid = rows[(size_t)U_indiv(rng)];
          for (int j = 0; j < p; ++j)
            mat_set_colmajor(mbuf.data(), n, rr, j,
                             mat_get_colmajor(dat.begin(), n, rid, j));
          rr++;
        }
      }

      wc84_fst_hs_ht_all_loci_ptr_ld(
        mbuf.data(), rr, n, p,
        pop_col, missing_code, base,
        fst_tmp.data(), hs_tmp.data(), ht_tmp.data()
      );

      const size_t row0 = (size_t)b * (size_t)L;
      for (int ell = 0; ell < L; ++ell) hs_boot_buf[row0 + (size_t)ell] = hs_tmp[(size_t)ell];
      hs_overall_buf[(size_t)b] = mean_across_loci_ptr(hs_tmp.data(), L);
    }
  }

  NumericMatrix hs_boot(B, L);
  NumericVector hs_overall_boot(B, NA_REAL);
  for (int b = 0; b < B; ++b) {
    hs_overall_boot[b] = hs_overall_buf[(size_t)b];
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) hs_boot(b, ell) = hs_boot_buf[row0 + (size_t)ell];
  }

  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);
  colnames(hs_boot) = loc_names;

  return List::create(
    _["HS_boot"]         = hs_boot,
    _["HS_overall_boot"] = hs_overall_boot,
    _["locus_names"]     = loc_names,
    _["boot_type"]       = "individual",
    _["n_boot"]          = B
  );
}