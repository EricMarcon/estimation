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
          inputId = "sd", 
          label = "Lognormal Standard Deviation:", 
          value = 2, 
          min = 0
        )
      ),
      conditionalPanel(
        condition = "input.distribution == 'geom'",
        numericInput(
          inputId = "prob",
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
          inputId = "Falpha",
          label = "Fisher's alpha:", 
          value = 40, 
          min = 1
        )
      ),
      sliderInput(
        inputId = "richness",
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
        inputId = "samplesize",
        label = "Sample Size",
        min = 100,
        max = 10000,
        value = 5000
      ),
      actionButton(inputId = "run", label = "Run"),
      helpText("More options"),
      sliderInput(
        inputId = "nsimulations",
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
    
    ## Main panel ####
    mainPanel(
      tabsetPanel(
        tabPanel("Diversity", plotOutput("profile")),
        tabPanel("RMSE", plotOutput("rmse")),
        tabPanel("Coverage", plotOutput("C")),
        tabPanel("RAC", plotOutput("rac")),
        tabPanel(
          "Help",
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
            href = "https://ericmarcon.github.io/entropart/",
            tags$img(
              src = "entropart.png",
              title = "Made with the entropart package",
              style = "display: block; margin-left: auto; margin-right: auto;"
              )
            )
          )
        )
      )
    )
  )

# Server logic ####
server <- function(input, output) {
  ## Store the reactive values ####
  RAC <- reactiveVal()
  Profile <- reactiveVal()
  RMSE <- reactiveVal()
  Sims_C <- reactiveVal()

  ## Action launched by the button ####
  observeEvent(
    input$run,
    {
      # Create the community
      if (input$distribution %in% c("lnorm", "geom"))
        Size <- .Machine$integer.max
      if (input$distribution  %in% c("bstick"))
        Size <- .Machine$integer.max %/% 2
      if (input$distribution  %in% c("lseries"))
        Size <- .Machine$integer.max %/% 2^6
      Community <- rCommunity(
        1,
        size = Size,
        S = input$richness,
        Distribution = input$distribution,
        sd = input$sd,
        prob = input$prob,
        alpha = input$Falpha,
        CheckArguments = FALSE
        )
      RAC(Community)
      Ps <- as.ProbaVector(Community)
      # Real profile
      q.seq <- c(seq(0, .1, 0.025), .2, seq(.3, 2, 0.1))
      alpha <- input$alpha
        
      # Compute the profile
      withProgress({
        # Estimated profile
        Values <- vapply(
          q.seq, 
          function(q) 
            Diversity(
              Ps,
              q,
              Correction = input$estimator,
              CheckArguments = FALSE
            ),
          0)
        
        # Create a MetaCommunity made of simulated communities
        MCSim <- rCommunity(
          input$nsimulations,
          size = input$samplesize,
          NorP = Ps,
          sd = input$sd,
          prob = input$prob,
          CheckArguments = FALSE
        )
        
        # Sample coverages of simulated communities: real and estimated
        Sims_C_real <- apply(
          MCSim$Nsi, 
          2,
          function(community)
            sum(Ps[community > 0])
          )
        Sims_C_est <- apply(MCSim$Nsi, 2, Coverage, CheckArguments = FALSE)
        sims_c <- tibble(
          Actual = Sims_C_real,
          Estimated = Sims_C_est
          )
        # Save the values
        Sims_C(sims_c)
        
        # Prepare a matrix for simulation results
        Sims <- matrix(
          nrow = input$nsimulations,
          ncol = length(q.seq)
          )
        
        # Run the simulations
        for (i in 1:input$nsimulations) {
          # Parallelize. Do not allow more forks in PhyloApply()
          ProfileAsaList <- parallel::mclapply(
            q.seq, 
            function(q)
              Diversity(
                MCSim$Nsi[, i],
                q,
                Correction = input$estimator,
                CheckArguments = FALSE
              ),
            mc.allow.recursive = FALSE)
          Sims[i, ] <- simplify2array(ProfileAsaList)
          setProgress(i)
        }
        
        Means <- apply(Sims, 2, mean)
        Vars <- apply(Sims, 2, var)
        # Quantiles of simulations for each q
        EstEnvelope <- apply(
          Sims,
          2,
          stats::quantile,
          probs = c(alpha / 2, 1 - alpha / 2)
          )
        colnames(EstEnvelope) <- q.seq
        cprofile <- list(
          x = q.seq,
          y = Values,
          low = EstEnvelope[1, ],
          high = EstEnvelope[2, ],
          mid = Means,
          var = Vars
          )
        class(cprofile) <- "CommunityProfile"
      },
      min = 0,
      max = input$nsimulations,
      message = "Running Simulations"
      )
      # Save the profile
      Profile(cprofile)
      # Calculate RMSE
      rmse <- tibble(
        q = q.seq,
        RMSE = sqrt((cprofile$y - cprofile$mid) ^ 2 + cprofile$var) / cprofile$y
      )
      RMSE(rmse)
    })
    
    # Output ####
    output$profile <- renderPlot(
      {
        if (inherits(Profile(), what = "CommunityProfile"))
          autoplot(Profile()) +
          labs(
            title = "Estimated vs Actual Diversity Profile",
            caption = paste(
              "Actual diversity profile of the comunity (solid line) and estimated diversity (confidence envelope at the",
              input$alpha,
              "risk level).
              The dotted line is the average estimation accross",
              input$nsimulations,
              "simulations"
              )
            )
        }
      )
    output$rmse <- renderPlot({
      if (inherits(Profile(), what = "CommunityProfile")) {
        rmse <- RMSE()
        ggplot(rmse) +
          geom_line(aes(x = q, y = RMSE)) +
          labs(
            title = "Normalized Root-Mean-Square Deviation",
                 caption = "RMSE normalized by the average value (accross simulations) of the estimator."
            )
        }
      })
    output$C <- renderPlot({
      if (inherits(Profile(), what = "CommunityProfile")) {
        sims_c <- Sims_C()
        ggplot(sims_c) + 
          geom_point(aes(x = Actual, y = Estimated)) + 
          geom_abline(col = "red") + 
          geom_hline(yintercept = mean(sims_c$Estimated), lty = 2) + 
          geom_vline(xintercept = mean(sims_c$Actual), lty = 2) +
          labs(
            title = "Sample Coverage",
            caption = "Estimated vs actual sample coverage of simulated samples.
            The red line represents equality.
            The dotted lines show the average estimated and actual values."
          )
        }
    })
    output$rac <- renderPlot({
      if (inherits(Profile(), what = "CommunityProfile")) {
        autoplot(RAC(), Distribution = input$distribution) +
          labs(
            title = "Rank-Abundance Curve (Whittaker Plot) of the simulated community.",
            caption = "Species are ranked from the most abundant to the rarest. 
            The abundance axis is in log-scale.
            The red curve is the best fit of the theoretical distribution."
            )
        }
    })
}

# Prepare the application ####

# Does the app run locally?
is_local <- (Sys.getenv('SHINY_PORT') == "")

# Install necessary packages ####
InstallPackages <- function(Packages) {
  sapply(Packages, function(Package)
    if (!Package %in% installed.packages()[, 1]) {
      install.packages(Package)
    })
}

# Necessary packages (not run on shyniapps.io)
if (is_local)
  InstallPackages(c("shiny", "tidyverse", "entropart", "parallel"))

# Load packages
library("shiny")
library("tidyverse")
library("entropart")


# Run the application ####
shinyApp(ui = ui, server = server)
