#' Splash Screen UI
#'
#' Fully isolated landing page. All styles are scoped under .spg-splash-root
#' and injected/removed dynamically via JavaScript when Launch is clicked.
#' No global body/html rules — zero bleed into app_ui().
#'
#' @return A \code{shiny::tagList} object.
#' @export
splash_ui <- function() {

  shiny::addResourcePath(
    "spg_www",
    system.file("app/www", package = "shinypopgen")
  )

  logos <- spg_logo_uris()

  # ── Helpers ───────────────────────────────────────────────────────────────
  logo_pill <- function(href, img_src, alt, fallback_text,
                        fallback_color = "#FFFFFF") {
    content <- if (!is.null(img_src))
      shiny::tags$img(
        src   = img_src, alt = alt,
        style = "max-height:36px;max-width:140px;object-fit:contain;
                 filter:brightness(1.15);"
      )
    else
      shiny::tags$span(
        style = paste0("font-size:13px;font-weight:700;letter-spacing:.5px;
                        color:", fallback_color, ";"),
        fallback_text
      )
    shiny::tags$a(
      href   = href, target = "_blank",
      style  = "text-decoration:none;",
      shiny::div(class = "spg-logo-pill", content)
    )
  }

  stat_card <- function(icon_name, value, label) {
    shiny::div(
      class = "spg-stat-card",
      shiny::icon(icon_name,
        style = "font-size:1.1rem;color:rgba(79,195,247,0.35);
                 float:right;margin-top:-2px;"
      ),
      shiny::div(class = "spg-stat-num",   value),
      shiny::div(class = "spg-stat-lbl",   label)
    )
  }

  mod_chip <- function(icon_name, label) {
    shiny::div(
      class = "spg-chip",
      shiny::icon(icon_name), label
    )
  }

  # ── CSS — ALL rules scoped to .spg-splash-root ────────────────────────────
  # Nothing touches html, body, or any global selector.
  # The <style> tag itself carries id="spg-splash-style" so JS can remove it.
  splash_css <- shiny::tags$style(
    id = "spg-splash-style",
    shiny::HTML("
/* ── Wrapper fills the uiOutput container ── */
.spg-splash-root {
  display: grid;
  grid-template-rows: auto 1fr auto;
  min-height: 100vh;
  width: 100%;
  background: #0b1a3d;
  font-family: 'Helvetica Neue','Segoe UI',Arial,sans-serif;
  color: #fff;
  box-sizing: border-box;
  overflow-x: hidden;
}
.spg-splash-root *, .spg-splash-root *::before, .spg-splash-root *::after {
  box-sizing: border-box;
}

/* top bar */
.spg-top-bar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 20px 48px 0;
}
.spg-logo-mark   { display:flex;align-items:center;gap:10px; }
.spg-logo-hex    {
  width:38px;height:38px;border-radius:8px;
  background:rgba(79,195,247,0.12);
  border:1.5px solid rgba(79,195,247,0.35);
  display:flex;align-items:center;justify-content:center;
  font-size:1.1rem;color:#4fc3f7;
}
.spg-appname { font-size:17px;font-weight:600;color:#e8f4fd;letter-spacing:.3px; }
.spg-doc-btn {
  display:flex!important;align-items:center!important;gap:7px!important;
  background:rgba(0,229,200,0.10)!important;
  border:1.5px solid rgba(0,229,200,0.38)!important;
  border-radius:6px!important;
  color:#00e5c8!important;
  font-size:13px!important;font-weight:600!important;
  padding:7px 16px!important;
  text-decoration:none!important;
  transition:background .18s;
}
.spg-doc-btn:hover { background:rgba(0,229,200,0.20)!important; }

/* main 2-col */
.spg-main {
  display:grid;
  grid-template-columns:1fr 1fr;
  padding:36px 48px 28px;
  align-items:center;
  gap:0;
}
@media(max-width:860px){
  .spg-main { grid-template-columns:1fr; }
  .spg-hero-right { padding-left:0!important;margin-top:28px; }
}

/* hero left */
.spg-badge {
  display:inline-flex;align-items:center;gap:6px;
  background:rgba(79,195,247,0.10);
  border:1px solid rgba(79,195,247,0.26);
  border-radius:20px;
  font-size:11px;font-weight:600;color:#4fc3f7;
  padding:4px 13px;letter-spacing:.8px;text-transform:uppercase;
  margin-bottom:16px;
}
.spg-hero-title {
  font-size:2.6rem;font-weight:700;color:#fff;
  line-height:1.15;margin-bottom:4px;
}
.spg-hero-title .spg-accent { color:#4fc3f7; }
.spg-hero-version {
  font-size:.95rem;color:rgba(79,195,247,.65);
  font-weight:500;margin-bottom:18px;letter-spacing:.5px;
}
.spg-hero-desc {
  font-size:1rem;color:#b8d4ef;line-height:1.65;
  margin-bottom:24px;max-width:420px;
}
.spg-authors-block {
  background:rgba(255,255,255,0.05);
  border:1px solid rgba(255,255,255,0.09);
  border-radius:8px;padding:14px 18px;
  margin-bottom:30px;max-width:420px;
}
.spg-authors-block .arow {
  display:flex;align-items:baseline;gap:8px;
  font-size:13px;color:#8ab4d4;line-height:1.9;
}
.spg-authors-block .arow strong { color:#c8dff0;font-weight:600; }
.spg-authors-block .ainst {
  font-size:12px;color:rgba(138,180,212,.68);
  margin-top:6px;padding-top:8px;
  border-top:1px solid rgba(255,255,255,0.07);
}

/* launch button */
.spg-launch-btn,
.spg-launch-btn.btn,
.spg-launch-btn.btn-default {
  display:inline-flex!important;align-items:center!important;gap:10px!important;
  background:#00e5c8!important;
  border:none!important;border-radius:8px!important;
  color:#0b1a3d!important;font-size:1rem!important;font-weight:700!important;
  padding:13px 36px!important;
  box-shadow:0 0 28px rgba(0,229,200,0.16);
  transition:background .18s,transform .12s;
}
.spg-launch-btn:hover,
.spg-launch-btn.btn:hover {
  background:#00ffdc!important;
  transform:translateY(-2px);
  color:#0b1a3d!important;
}

/* hero right */
.spg-hero-right {
  display:flex;flex-direction:column;align-items:center;
  gap:22px;padding-left:28px;
}
.spg-stats-grid {
  display:grid;grid-template-columns:1fr 1fr;
  gap:12px;width:100%;max-width:340px;
}
.spg-stat-card {
  background:rgba(255,255,255,0.05);
  border:1px solid rgba(255,255,255,0.10);
  border-radius:10px;padding:16px 18px;
}
.spg-stat-num { font-size:1.55rem;font-weight:700;color:#4fc3f7;margin-bottom:3px; }
.spg-stat-lbl {
  font-size:11px;color:#8ab4d4;font-weight:600;
  text-transform:uppercase;letter-spacing:.7px;
}
.spg-module-list { width:100%;max-width:340px; }
.spg-module-list .spg-mtitle {
  font-size:11px;font-weight:600;color:#5a8aaa;
  text-transform:uppercase;letter-spacing:.8px;margin-bottom:10px;
}
.spg-chips { display:flex;flex-wrap:wrap;gap:7px; }
.spg-chip {
  background:rgba(255,255,255,0.07);
  border:1px solid rgba(255,255,255,0.12);
  border-radius:5px;font-size:12px;color:#b8d4ef;
  padding:5px 11px;display:flex;align-items:center;gap:5px;
}
.spg-chip .fa { font-size:.78rem;color:rgba(79,195,247,0.55); }

/* logos bar */
.spg-logos-bar {
  background:rgba(0,0,0,0.22);
  border-top:1px solid rgba(255,255,255,0.07);
  padding:18px 48px;
  display:flex;align-items:center;
  justify-content:space-between;gap:16px;flex-wrap:wrap;
}
.spg-logos-label {
  font-size:11px;color:rgba(138,180,212,.52);
  font-weight:600;letter-spacing:.9px;text-transform:uppercase;
  white-space:nowrap;flex-shrink:0;
}
.spg-logo-items { display:flex;align-items:center;gap:20px;flex-wrap:wrap; }
.spg-logo-pill {
  background:rgba(255,255,255,0.09);
  border:1px solid rgba(255,255,255,0.14);
  border-radius:7px;padding:8px 18px;
  display:flex;align-items:center;justify-content:center;
  min-height:48px;min-width:80px;transition:background .18s;
}
.spg-logo-pill:hover { background:rgba(255,255,255,0.17); }
    ")
  )

  # ── JS: on Launch click, remove splash CSS + the root div ────────────────
  # Shiny's observeEvent in app_server handles the actual ui swap.
  # This JS only cleans up the style tag so no splash rule bleeds into app_ui.
  cleanup_js <- shiny::tags$script(shiny::HTML("
$(document).on('click', '#btn_launch', function() {
  // Remove the scoped splash stylesheet immediately
  var s = document.getElementById('spg-splash-style');
  if (s) s.parentNode.removeChild(s);
});
  "))

  # ── Markup ───────────────────────────────────────────────────────────────
  shiny::tagList(
    splash_css,
    cleanup_js,

    shiny::div(
      class = "spg-splash-root",

      # TOP BAR
      shiny::div(
        class = "spg-top-bar",
        shiny::div(
          class = "spg-logo-mark",
          shiny::div(class = "spg-logo-hex", shiny::icon("dna")),
          shiny::span(class = "spg-appname", "ShinyPopGen")
        ),
        shiny::tags$a(
          class  = "spg-doc-btn",
          href   = "https://forge.ird.fr/intertryp/shiny_pop_gen",
          target = "_blank",
          shiny::icon("file-alt"), " Documentation"
        )
      ),

      # MAIN
      shiny::div(
        class = "spg-main",

        # Left
        shiny::div(
          shiny::div(class = "spg-badge", shiny::icon("flask"),
                     " Population genetics"),
          shiny::div(
            class = "spg-hero-title",
            "ShinyPopGen", shiny::tags$br(),
            shiny::tags$span(class = "spg-accent", "V1")
          ),
          shiny::div(class = "spg-hero-version",
                     "SPG \u00b7 Interactive Analysis Suite"),
          shiny::div(
            class = "spg-hero-desc",
            "A versatile, user-friendly and multi-OS application to analyse
             population genetic data \u2014 import, explore, and run descriptive
             statistics in a few clicks."
          ),
          shiny::div(
            class = "spg-authors-block",
            shiny::div(class = "arow",
              shiny::tags$strong("Programming"),
              " Vincent Manzanilla & Naffiou Kaderi"),
            shiny::div(class = "arow",
              shiny::tags$strong("Conception"),
              " Thierry de Mee\u00fbs"),
            shiny::div(class = "ainst",
              "Intertryp \u00b7 Univ. Montpellier \u00b7 Cirad \u00b7
               IRD \u00b7 Montpellier, France")
          ),
          shiny::actionButton(
            inputId = "btn_launch",
            label   = shiny::tagList(shiny::icon("rocket"), " Launch application"),
            class   = "spg-launch-btn"
          )
        ),

        # Right
        shiny::div(
          class = "spg-hero-right",
          shiny::div(
            class = "spg-stats-grid",
            stat_card("th-large",        "10",     "Analysis modules"),
            stat_card("upload",          "500 MB", "Upload limit"),
            stat_card("project-diagram", "WC84",   "F-statistics"),
            stat_card("desktop",         "Multi",  "OS compatible")
          ),
          shiny::div(
            class = "spg-module-list",
            shiny::div(class = "spg-mtitle", "Available modules"),
            shiny::div(
              class = "spg-chips",
              mod_chip("upload",         "Data import"),
              mod_chip("chart-pie",      "Allele freq."),
              mod_chip("table",          "General stats"),
              mod_chip("flask",          "Local panmixia"),
              mod_chip("globe",          "Global panmixia"),
              mod_chip("sitemap",        "Subdivision"),
              mod_chip("chart-line",     "Diversities"),
              mod_chip("link",           "LD"),
              mod_chip("circle-notch",   "Null alleles"),
              mod_chip("map-marker-alt", "IBD")
            )
          )
        )
      ),

      # LOGOS BAR
      shiny::div(
        class = "spg-logos-bar",
        shiny::span(class = "spg-logos-label", "Partners"),
        shiny::div(
          class = "spg-logo-items",
          logo_pill("https://umr-intertryp.cirad.fr/en",
                    logos$intertryp, "Intertryp", "INTERTRYP"),
          logo_pill("https://www.ird.fr/en",
                    logos$ird, "IRD", "IRD", "#e8c44a"),
          logo_pill("https://www.umontpellier.fr/en/",
                    logos[["um"]], "Université de Montpellier",
                    "UNIVERSIT\u00c9 DE MONTPELLIER", "#e05050"),
          logo_pill("https://www.cirad.fr/en",
                    logos$cirad, "CIRAD", "CIRAD", "#7cc576")
        )
      )
    )
  )
}