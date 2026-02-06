# UI ----
ui <- fluidPage(
  # Application title
  titlePanel("Diversity Estimation"),

  sidebarLayout(
    ## Sidebar ----
    sidebarPanel(
      selectInput(
        inputId = "distribution",
        label = "Distribution:",
        choices = c(
          "Lognormal" = "lnorm",
          "Geometric" = "geom",
          "Log-series" = "lseries",
          "Broken Stick" = "bstick"
        )
      ),
      conditionalPanel(
        condition = "input.distribution == 'lnorm'",
        numericInput(
          inputId = "sd_lnorm",
          label = "Lognormal Standard Deviation:",
          value = 2,
          min = 0
        )
      ),
      conditionalPanel(
        condition = "input.distribution == 'geom'",
        numericInput(
          inputId = "prob_geom",
          label = "Ressource Share:",
          value = 0.01,
          min = 0,
          max = 0.5,
          step = 0.01
        )
      ),
      conditionalPanel(
        condition = "input.distribution == 'lseries'",
        numericInput(
          inputId = "fisher_alpha",
          label = "Fisher's alpha:",
          value = 40,
          min = 1
        )
      ),
      sliderInput(
        inputId = "species_number",
        label = "Number of Species in the Community",
        min = 2,
        max = 1000,
        value = 300
      ),
      selectInput(
        inputId = "estimator",
        label = "Estimator:",
        choices = c(
          "Chao-Jost" = "ChaoJost",
          "Chao-Shen" = "ChaoShen",
          "Generalized Coverage" = "GenCov",
          "Unveiled Chao1" = "UnveilC",
          "Unveiled Improved Chao1" = "UnveiliC",
          "Unveiled Jackknife" = "UnveilJ"
        ),
        selected = "UnveilJ"
      ),
      sliderInput(
        inputId = "sample_size",
        label = "Sample Size",
        min = 100,
        max = 10000,
        value = 5000
      ),
      actionButton(inputId = "run", label = "Run"),
      helpText("More options"),
      sliderInput(
        inputId = "n_simulations",
        label = "Number of Simulations",
        min = 2,
        max = 1000,
        value = 10
      ),
      numericInput(
        inputId = "alpha",
        label = "Risk level:",
        value = 0.05,
        min = 0.001,
        max = 0.5,
        step = 0.01
      ),
    ),

    ## Main panel ----
    mainPanel(
      tabsetPanel(
        tabPanel(title = "Diversity", plotOutput("profile")),
        tabPanel(title = "RMSE", plotOutput("rmse")),
        tabPanel(title = "Coverage", plotOutput("C")),
        tabPanel(title = "RAC", plotOutput("rac")),
        tabPanel(
          title = "Help",
          p(),
          p(
            "This Shiny app is designed to test the efficiency of various diversity estimators applied to undersampled communities",
            a(href = "https://hal-agroparistech.archives-ouvertes.fr/hal-01212435v2/document", "(Marcon, 2015)"),
            "."
            ),
          p(
            "It generates a community from a well-known distribution, simulates sampling and estimates diversity from the samples.",
            "The estimated diversity is compared to the actual one."
            ),
          h2("Generate a community"),
          p(
            "Choose the distribution and its parameters, including the number of species.",
            "A community of the largest size allowed by R is created (only one billion individuals for the broken-stick distribution)."
            ),
          p(
            "Its Rank-Abundance Curve, aka a Whittaker plot, is displayed in the RAC tab of the app."
            ),
          h2("Sample it"),
          p(
            "Choose the sample size and the number of simulations.",
            "10 simulations are enough for a quick view of the estimation.",
            "Run 1000 simulations for accurate statistics."
            ),
          h2("Estimate"),
          p("Choose the estimator of diversity."),
          p(
            "It is applied to each sample to obtain a distribution of estimations."
          ),
          p(
            "The confidence envelope of diversity profiles compared to the diversity of the actual community is in the Diversity tab."
          ),
          p(
            "The Root-Mean-Square Deviation normalized by the average estimated diversity is in the RMSE tab."
          ),
          p(
            "The estimated and actual sample coverage of each simulated sample are in the Coverage tab."
          ),
          tags$a(
            href = "https://ericmarcon.github.io/divent/",
            tags$img(
              src = "divent.png",
              title = "Made with the divent package",
              style = "display: block; margin-left: auto; margin-right: auto;"
              )
            )
          )
        )
      )
    )
  )

# Server logic ----
server <- function(input, output) {
  ## Store the reactive values ----
  rac_reactive <- reactiveVal()
  profile_reactive <- reactiveVal()
  rmse_reactive <- reactiveVal()
  coverages_compare_reactive <- reactiveVal()

  ## Action launched by the button ----
  observeEvent(
    eventExpr = input$run,
    handlerExpr = {
      # Create the community
      if (input$distribution %in% c("lnorm", "geom"))
        the_size <- .Machine$integer.max
      if (input$distribution %in% c("bstick"))
        the_size <- .Machine$integer.max %/% 2
      if (input$distribution %in% c("lseries"))
        the_size <- .Machine$integer.max %/% 2^6
      the_community <- rcommunity(
        n = 1,
        size = the_size,
        species_number = input$species_number,
        distribution = input$distribution,
        sd_lnorm = input$sd_lnorm,
        prob_geom = input$prob_geom,
        fisher_alpha = input$fisher_alpha,
        check_arguments = FALSE
        )
      rac_reactive(the_community)
      prob <- as.numeric(as_probabilities(the_community))
      # Real profile
      orders <- c(seq(0, .1, 0.025), .2, seq(.3, 2, 0.1))
      alpha <- input$alpha

      # Compute the profile
      withProgress({
        # Estimated profile
        the_diversity <- vapply(
          orders,
          function(q) {
            div_hill(
              prob,
              q = q,
              as_numeric = TRUE,
              check_arguments = FALSE
            )
          },
          FUN.VALUE = 0
        )

        # Create a MetaCommunity made of simulated communities
        the_metacommunity <- rcommunity(
          n = input$n_simulations ,
          size = input$sample_size,
          prob = prob,
          check_arguments = FALSE
        )

        # Sample coverage of simulated communities: real and estimated
        coverages_real <- apply(
          as.matrix(the_metacommunity),
          MARGIN = 1,
          function(community)
            sum(prob[community > 0])
          )
        coverages_estimated <- coverage(
          the_metacommunity,
          as_numeric = TRUE,
          check_arguments = FALSE
        )
        coverages_compare <- tibble(
          Actual = coverages_real,
          Estimated = coverages_estimated
          )
        # Save the values
        coverages_compare_reactive(coverages_compare)

        the_metacommunity |>
          profile_hill(orders = orders) |>
          pivot_wider(
            names_from = site,
            values_from = diversity
          ) |>
          select(-estimator, -order) ->
          the_metacommunity_hill

        the_metacommunity_hill_mean <- rowMeans(the_metacommunity_hill)
        the_metacommunity_hill_var <- apply(
          the_metacommunity_hill,
          MARGIN = 1,
          FUN = var
        )
        # Quantiles of simulations for each q
        the_metacommunity_hill_envelope <- apply(
          the_metacommunity_hill,
          MARGIN = 1,
          FUN = stats::quantile,
          probs = c(alpha / 2, 1 - alpha / 2)
          )
        colnames(the_metacommunity_hill_envelope) <- orders

        # Make a profile object to plot it
        the_profile <- tibble(
          site = "simulations",
          order = orders,
          diversity = the_diversity,
          inf = the_metacommunity_hill_envelope[1, ],
          sup = the_metacommunity_hill_envelope[2, ],
          # Extra fields
          mean = the_metacommunity_hill_mean,
          var = the_metacommunity_hill_var
        )
        class(the_profile) <- c(
          "profile",
          class(the_profile)
        )
      },
      min = 0,
      max = input$n_simulations,
      message = "Running Simulations"
      )
      # Save the profile
      profile_reactive(the_profile)
      # Calculate RMSE
      rmse <- tibble(
        q = orders,
        RMSE = sqrt((the_profile$diversity - the_profile$mean) ^ 2 + the_profile$var) / the_profile$diversity
      )
      rmse_reactive(rmse)
    })

    # Output ----
    output$profile <- renderPlot(
      {
        if (inherits(profile_reactive(), what = "profile"))
          autoplot(profile_reactive()) +
          labs(
            title = "Estimated vs Actual Diversity Profile",
            caption = paste(
              "Actual diversity profile of the comunity (solid line) and estimated diversity (confidence envelope at the",
              input$alpha,
              "risk level).
              The dotted line is the average estimation accross",
              input$n_simulations ,
              "simulations"
              )
            )
        }
      )
    output$rmse <- renderPlot({
      if (inherits(profile_reactive(), what = "profile")) {
        rmse <- rmse_reactive()
        ggplot(rmse) +
          geom_line(aes(x = q, y = RMSE)) +
          labs(
            title = "Normalized Root-Mean-Square Deviation",
                 caption = "RMSE normalized by the average value (accross simulations) of the estimator."
            )
        }
      })
    output$C <- renderPlot({
      if (inherits(profile_reactive(), what = "profile")) {
        coverages_compare <- coverages_compare_reactive()
        ggplot(coverages_compare) +
          geom_point(aes(x = Actual, y = Estimated)) +
          geom_abline(col = "red") +
          geom_hline(yintercept = mean(coverages_compare$Estimated), lty = 2) +
          geom_vline(xintercept = mean(coverages_compare$Actual), lty = 2) +
          labs(
            title = "Sample Coverage",
            caption = "Estimated vs actual sample coverage of simulated samples.
            The red line represents equality.
            The dotted lines show the average estimated and actual values."
          )
        }
    })
    output$rac <- renderPlot({
      if (inherits(profile_reactive(), what = "profile")) {
        autoplot(rac_reactive(), Distribution = input$distribution) +
          labs(
            title = "Rank-Abundance Curve (Whittaker Plot) of the simulated community.",
            caption = "Species are ranked from the most abundant to the rarest.
            The abundance axis is in log-scale.
            The red curve is the best fit of the theoretical distribution."
            )
        }
    })
}

# Prepare the application ----

# Does the app run locally?
is_local <- (Sys.getenv('SHINY_PORT') == "")

# Install necessary packages ----
install_packages <- function(packages) {
  invisible(
    sapply(
      packages,
      FUN = function(package) {
        if (!package %in% installed.packages()[, 1]) {
          install.packages(
            package,
            repos = "https://cran.rstudio.com/",
            quiet = TRUE
          )
        }
      }
    )
  )
}

# Necessary packages (not run on shyniapps.io)
if (is_local)
  install_packages(c("shiny", "tidyverse", "divent", "parallel"))

# Load packages
library("shiny")
library("tidyverse")
library("divent")

# Run the application ----
shinyApp(ui = ui, server = server)
