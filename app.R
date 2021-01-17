# UI ####
ui <- fluidPage(# Application title
    titlePanel("Diversity Estimation"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput(
                "distribution", 
                "Distribution:",
                c("Lognormal" = "lnorm",
                  "Geometric" = "geom")
                ),
            conditionalPanel(
                condition = "input.distribution == 'lnorm'",
                numericInput("sd", "Lognormal Standard Deviation:", 2, min = 0)
                ),
            conditionalPanel(
                condition = "input.distribution == 'geom'",
                numericInput("prob ", "Ressource Share:", 0.1, min = 0, max = 1, step=0.1)
                ),
            sliderInput(
                "richness",
                "Number of Species in the Community",
                min = 2,
                max = 1000,
                value = 300
                ),
            selectInput(
                "estimator", 
                "Estimator:",
                c("Unveiled Jackknife" = "UnveilJ",
                  "Chao-Jost" = "ChaoJost")
            ),
            sliderInput(
                "samplesize",
                "Sample Size",
                min = 100,
                max = 10000,
                value = 5000
            ),
            actionButton("run", "Run"),
            helpText("More options"),
            sliderInput(
                "nsimulations",
                "Number of Simulations",
                min = 2,
                max = 1000,
                value = 10
            ),
            numericInput("alpha", "Risk level:", 0.05, min = 0.001, max = 0.5, step=0.01),
            # End of input
        ),
        
        # Show plots in the main panel
        mainPanel(
            tabsetPanel(
                tabPanel("Diversity", plotOutput("profile")), 
                tabPanel("RMSE", plotOutput("rmse")), 
                tabPanel("RAC", plotOutput("rac"))
            )
        )
    ))

# Server logic ####
server <- function(input, output) {
    
    # Store the values to plot
    RAC <- reactiveVal()
    Profile <- reactiveVal()
    RMSE <- reactiveVal()
    
    # Action launched by the button
    observeEvent(input$run, {
        # Create the community
        Community <- rCommunity(1,                               
                                size=1E6,
                                S=input$richness,
                                Distribution=input$distribution,
                                CheckArguments = FALSE
                                )
        RAC(Community)
        Ps <- as.ProbaVector(Community)
        # Real profile
        q.seq <- c(seq(0, .1, 0.025), .15, seq(.2, 1, 0.1))
        Alpha <- input$alpha
        
        # Compute the profile
        withProgress({
            # Estimated profile
            Values <- vapply(q.seq, function(q) Diversity(Ps, q, CheckArguments = FALSE), 0)
            # Create a MetaCommunity made of simulated communities
            MCSim <- rCommunity(input$nsimulations, size=input$samplesize, NorP=Ps, CheckArguments = FALSE)
            Sims <- matrix(nrow=input$nsimulations, ncol=length(q.seq))
            for (i in 1:input$nsimulations) {
                # Parallelize. Do not allow more forks in PhyloApply()
                ProfileAsaList <- parallel::mclapply(q.seq, function(q) Diversity(MCSim$Nsi[, i], q, CheckArguments=FALSE), mc.allow.recursive=FALSE)
                Sims[i, ] <- simplify2array(ProfileAsaList)
                setProgress(i)
            }
            Means <- apply(Sims, 2, mean)
            Vars <- apply(Sims, 2, var)
            # Quantiles of simulations for each q
            EstEnvelope <- apply(Sims, 2, stats::quantile, probs = c(Alpha/2, 1-Alpha/2))
            colnames(EstEnvelope) <- q.seq
            cprofile <- list(x=q.seq,
                            y=Values,
                            low=EstEnvelope[1, ],
                            high=EstEnvelope[2, ],
                            mid=Means,
                            var=Vars)
            class(cprofile) <- "CommunityProfile"
            }, 
            min = 0,
            max = input$nsimulations,
            message = "Running Simulations"
        )
        # Save the profile
        Profile(cprofile)
        # Calculate RMSE
        rmse <- tibble(q=q.seq, RMSE=sqrt((cprofile$y-cprofile$mid)^2 + cprofile$var)/cprofile$y)
        RMSE(rmse)
    })
    
    # Output ####
    output$profile <- renderPlot({
        if (inherits(Profile(), what="CommunityProfile"))
            autoplot(Profile())
    })
    output$rmse <- renderPlot({
        if (inherits(Profile(), what="CommunityProfile")) {
            rmse <- RMSE()
            ggplot(rmse) +
                geom_line(aes(x=q, y=RMSE))
        }
    })
    output$rac <- renderPlot({
        if (inherits(Profile(), what="CommunityProfile")) {
            autoplot(RAC(), Distribution=input$distribution)
        }
    })
}

# Prepare the application ####

# Does the app run locally?
is_local <- (Sys.getenv('SHINY_PORT') == "")

# Install necessary packages ####
InstallPackages <- function(Packages) {
    sapply(Packages, function(Package) 
        if (!Package %in% installed.packages()[, 1]) {install.packages(Package)})
}

# Necessary packages (not run on shyniapps.io)
if(is_local) InstallPackages(c("shiny", "tidyverse", "entropart", "parallel"))
    
# Load packages
library("shiny")
library("tidyverse")
library("entropart")


# Run the application ####
shinyApp(ui = ui, server = server)
