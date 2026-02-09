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
        inputId = "n_species",
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
          "Naive" = "naive",
          "Unveiled Chao1" = "UnveilC",
          "Unveiled Improved Chao1" = "UnveiliC",
          "Unveiled Jackknife" = "UnveilJ",
          "Zhang-Grabchak" = "ZhangGrabchak"
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
        tabPanel(title = "Diversity", plotOutput("div_profile")),
        tabPanel(title = "RMSE", plotOutput("rmse")),
        tabPanel(title = "Coverage", plotOutput("coverage")),
        tabPanel(title = "RAC", plotOutput("rac")),
        tabPanel(
          title = "Help",
          p(),
          p(
            "This Shiny app is designed to test the efficiency of various",
            "diversity estimators applied to undersampled communities",
            a(
              href = "https://hal-agroparistech.archives-ouvertes.fr/hal-01212435v2/document",
              "(Marcon, 2015)"
            ),
            "."
          ),
          p(
            "It generates a community from a well-known distribution,",
            "simulates sampling and estimates diversity from the samples.",
            "The estimated diversity is compared to the actual one."
          ),
          h2("Generate a community"),
          p(
            "Choose the distribution and its parameters,",
            "including the number of species.",
            "A community of the largest size allowed by R is created,",
            "except for log-series whose size depends on Fisher's alpha",
            "and the number of species."
          ),
          p(
            "Its Rank-Abundance Curve, aka a Whittaker plot,",
            "is displayed in the RAC tab of the app."
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
            "The confidence envelope of diversity profiles compared to",
            "the diversity of the actual community is in the Diversity tab."
          ),
          p(
            "The Root-Mean-Square Deviation normalized by the average",
            "estimated diversity is in the RMSE tab."
          ),
          p(
            "The estimated and actual sample coverage of each",
            "simulated sample are in the Coverage tab."
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
  the_community.rV <- reactiveVal()
  the_profile.rV <- reactiveVal()
  rmse.rV <- reactiveVal()
  c_compare.rV <- reactiveVal()
  # Fixed parameter: orders of diversity for the profile
  orders <- c(seq(0, .1, 0.025), .2, seq(.3, 2, 0.1))

  ## Action launched by the button ----
  observeEvent(
    eventExpr = input$run,
    handlerExpr = {
      # Create the community
      if (input$distribution %in% c("lnorm", "geom")) {
        the_size <- .Machine$integer.max
      }
      if (input$distribution %in% c("bstick")) {
        the_size <- .Machine$integer.max %/% 2
      }
      if (input$distribution %in% c("lseries")) {
        # Size depends on Fisher's alpha and richness
        the_size <- round(
          exp(
            (input$n_species + input$fisher_alpha * log(input$fisher_alpha)) /
              input$fisher_alpha
          ) - input$fisher_alpha
        )
      }
      # Draw a very large community
      the_community <- rcommunity(
        n = 1,
        size = the_size,
        species_number = input$n_species,
        distribution = input$distribution,
        sd_lnorm = input$sd_lnorm,
        prob_geom = input$prob_geom,
        fisher_alpha = input$fisher_alpha,
        check_arguments = FALSE
      )
      the_community.rV(the_community)
      # Store the probabilites of species
      prob <- as.numeric(as_probabilities(the_community))

      # Compute the profile
      withProgress(
        {
          # Diversity profile of the whole community
          the_diversity <- profile_hill(
            prob,
            orders = orders,
            as_numeric = TRUE
          )

          # Create a metacommunity made of simulated communities
          the_metacommunity <- rcommunity(
            n = input$n_simulations,
            size = input$sample_size,
            prob = prob,
            check_arguments = FALSE
          )
          # Store it in a matrix
          the_communities <- as.matrix(the_metacommunity)

          # Sample coverage of simulated communities: real and estimated
          c_real <- apply(
            the_communities,
            MARGIN = 1,
            FUN = function(community) {
              sum(prob[community > 0])
            }
          )
          c_estimated <- coverage(
            the_metacommunity,
            as_numeric = TRUE,
            check_arguments = FALSE
          )
          c_compare <- tibble(
            Actual = c_real,
            Estimated = c_estimated
          )
          # Save the values
          c_compare.rV(c_compare)

          # Prepare a matrix for profiles of simulated communities
          profile.matrix <- matrix(
            nrow = input$n_simulations,
            ncol = length(orders)
          )

          # Run the simulations separately to have a progress bar
          # Do profile_hill(the_metacommunity) equivalent
          for (i in seq_len(input$n_simulations)) {
            profile.matrix[i, ] <- profile_hill(
              the_communities[i, ],
              orders = orders,
              estimator = input$estimator,
              as_numeric = TRUE,
              check_arguments = FALSE
            )
            setProgress(i)
          }

          profile_mean <- colMeans(profile.matrix)
          profile_var <- apply(
            profile.matrix,
            MARGIN = 2,
            FUN = var
          )
          # Quantiles of simulations for each order
          profile_envelope <- apply(
            profile.matrix,
            MARGIN = 2,
            FUN = stats::quantile,
            probs = c(input$alpha / 2, 1 - input$alpha / 2)
            )

          # Make a profile object to plot it
          the_profile <- tibble(
            site = "simulations",
            order = orders,
            diversity = the_diversity,
            inf = profile_envelope[1, ],
            sup = profile_envelope[2, ],
            # Extra fields
            mean = profile_mean,
            var = profile_var
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
      the_profile.rV(the_profile)
      # Calculate RMSE
      rmse <- tibble(
        q = orders,
        RMSE = sqrt((the_profile$mean - the_profile$diversity) ^ 2 +
          the_profile$var) / the_profile$diversity
      )
      rmse.rV(rmse)
    })

    # Output ----
    output$div_profile <- renderPlot({
      if (inherits(the_profile.rV(), what = "profile")) {
        autoplot(the_profile.rV()) +
          theme(legend.position = "none") +
          geom_line(aes(y = mean), color = "black", lty = 2) +
          labs(
            title = "Estimated vs Actual Diversity Profile",
            caption = stringr::str_wrap(
              paste(
                "Actual diversity profile of the comunity (red line).",
                "The dotted line is the average estimation across",
                input$n_simulations,
                "simulations with its confidence envelope at the",
                input$alpha,
                "risk level)."
              )
            )
          )
      }
    })
    output$rmse <- renderPlot({
      if (inherits(the_profile.rV(), what = "profile")) {
        rmse <- rmse.rV()
        ggplot(rmse) +
          geom_line(aes(x = q, y = RMSE)) +
          labs(
            title = "Normalized Root-Mean-Square Deviation",
            caption = stringr::str_wrap(
              "RMSE normalized by the average value
              (across simulations) of the estimator."
            )
          )
      }
    })
    output$coverage <- renderPlot({
      if (inherits(the_profile.rV(), what = "profile")) {
        coverages_compare <- c_compare.rV()
        ggplot(coverages_compare) +
          geom_point(aes(x = Actual, y = Estimated)) +
          geom_abline(col = "red") +
          geom_hline(yintercept = mean(coverages_compare$Estimated), lty = 2) +
          geom_vline(xintercept = mean(coverages_compare$Actual), lty = 2) +
          labs(
            title = "Sample Coverage",
            caption = stringr::str_wrap(
              "Estimated vs actual sample coverage of simulated samples.
              The red line represents equality.
              The dotted lines show the average estimated and actual values."
            )
          )
      }
    })
    output$rac <- renderPlot({
      if (inherits(the_profile.rV(), what = "profile")) {
        autoplot(the_community.rV(), Distribution = input$distribution) +
          labs(
            title = "Rank-Abundance Curve (Whittaker Plot) of the simulated community.",
            caption = stringr::str_wrap(
              "Species are ranked from the most abundant to the rarest.
              The abundance axis is in log-scale.
              The red curve is the best fit of the theoretical distribution."
            )
          )
      }
    })
}

# Prepare the application ----

## Does the app run locally? ----
is_local <- (Sys.getenv('SHINY_PORT') == "")

## Install necessary packages ----
install_packages <- function(packages) {
  packages |>
    # Which packages are not installed yet?
    setdiff(installed.packages()[, 1]) |>
    # Install them
    install.packages(
      repos = "https://cran.rstudio.com/",
      quiet = TRUE
    ) |>
    # Hide the result
    invisible()
}

# Necessary packages (not run on shyniapps.io)
if (is_local) {
  install_packages(
    c("shiny", "tibble", "ggplot2", "divent")
  )
}


# Load packages
library("shiny")
library("tibble")
library("ggplot2")
library("divent")

# Run the application ----
shinyApp(ui = ui, server = server)
