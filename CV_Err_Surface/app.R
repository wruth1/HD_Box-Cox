library(shiny)
library(glmnet)
# library(stringr)
# library(ggplot2)
# library(dplyr)
# library(gtable)
# library(plot3D)
# library(reshape2)
library(plotly)

source("Helper Functions/All Helper Scripts.R")


ui <- fluidPage(
    # Application title
    titlePanel("CV Error Surface"),
    
    fluidRow(column(5,
                    numericInput("p",
                                 "p",
                                 value = 50)),
             column(
                 5,
                 # selectInput(
                 #     "q.str",
                 #     h3("q"),
                 #     choices = list("Full" = "full",
                 #                    "Sqrt" = "sqrt"),
                 #     selected = 1
                 # )
                 numericInput(
                     "q",
                     "Number of Active Predictors",
                     value = 5
                 )
             )),
    
    fluidRow(column(5,
                    numericInput(
                        "sigma",
                        "SD of Y",
                        value = 1
                    )),
             column(5,
                    numericInput(
                        "delta",
                        "SD of X %*% beta",
                        value = 1
                    ))),
    
    fluidRow(column(5,
                    numericInput(
                        "gamma.0",
                        "True Gamma",
                        value = 0
                    )),
             
             column(
                 5,
                 selectInput(
                     "CV.Type",
                     "Cross Validation Type",
                     choices = list("Min" = "min",
                                    "1SE" = "1se"),
                     selected = 1
                 )
             )),
    fluidRow(column(5,
                    numericInput(
                        "n",
                        "Sample Size",
                        value = 100
                    )),
             column(5,
                    checkboxGroupInput("to.log",
                                  "Variables to log:",
                                  choices = list("lambda" = "lambda",
                                                 "CV Error" = "response"),
                                  selected = c("lambda", "response")))),
    fluidRow(column(1,
                    actionButton("Button", "Plot")),
             column(9, "Please allow about some time to generate the plot")),
    fluidRow(plotlyOutput("CV.Surface.Plot"))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    set.seed(52314432)
    
    n = 100     #Sample Size
    n.lambda = 200 # Target number of lambda candidates
    
    ### Amount to decrease likelihood by to construct CI
    CI.step.size = qchisq(0.95, 1) / 2
    
    #Step size for gamma candidates
    gamma.grid.step = 0.01
    #Maximum increase and decrease from gamma.0
    gamma.step.up = 1
    gamma.step.down = 1
    
    ### Construct folds ###
    n.folds = 10
    fold.size = n %/% n.folds
    n.leftover = n %% n.folds
    folds.raw = rep(1:n.folds, times = fold.size)
    leftovers = seq_len(n.leftover)
    folds.raw = c(folds.raw, leftovers)
    folds = sample(folds.raw)
    
    
    output$CV.Surface.Plot <- renderPlotly({
        if (input$Button == 0) {
            return()
        }
        
        isolate({
            
            p = input$p
            # q.str = input$q.str
            q = input$q
            sigma = input$sigma
            delta = input$delta
            gamma.0 = input$gamma.0
            CV.Type = input$CV.Type
            to.log = input$to.log
            

            #Smallest and largest gamma candidates
            gamma.min = gamma.0 - gamma.step.down
            gamma.max = gamma.0 + gamma.step.up
            #Candidate gamma values
            Gammas = seq(gamma.min, gamma.max, gamma.grid.step)

            #lambda min or 1se
            lambda.type = paste0("lambda.", CV.Type)

            ### Construct coefficient vector
            # q = ifelse(q.str == "sqrt", sqrt(p), p)
            # q = floor(q)
            # q = floor(sqrt(p))
            beta.size = delta / sqrt(q)
            beta = c(rep(beta.size, q),
                     rep(0, p - q))



            ### Generate data
            X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
            mu.Z.raw = X %*% beta
            mu.Z = mu.Z.raw + get.int(n, delta, sigma)
            Z = mu.Z + rnorm(n, 0, sigma)

            ### Transform Z so that gamma.0 is correct BC parameter
            Y = inv.BC(Z, gamma.0)

            ### Find all candidate lambda values
            all.lambdas.raw = sapply(Gammas, function(gamma) {
                this.Z = BC(Y, gamma)
                this.fit = glmnet(X, this.Z)
                this.lambdas = this.fit$lambda
                return(this.lambdas)
            })
            all.lambdas.fine = sort(unlist(all.lambdas.raw))

            ### Make the grid of lambda candidates coarser
            all.lambdas = coarsen.grid(n.lambda, all.lambdas.fine)


            all.CV.errors = sapply(Gammas, function(gamma) {
                this.Z = BC(Y, gamma)
                this.fit = cv.glmnet(X,
                                     this.Z,
                                     lambda = all.lambdas,
                                     foldid = folds)
                this.errs = this.fit$cvm
                return(this.errs)
            })

            ### Pre-compute booleans for which variables to log
            log.lambda = "lambda" %in% to.log
            log.response = "response" %in% to.log
            
            if(!log.lambda & !log.response){
                plot_ly(
                    x = Gammas,
                    y = all.lambdas,
                    z = all.CV.errors
                ) %>%
                    add_surface() %>%
                    layout(scene = list(
                        xaxis = list(title = "gamma"),
                        yaxis = list(title = "lambda"),
                        zaxis = list(title = "CV-Error")
                    ))
            } else if(!log.lambda & log.response){
                plot_ly(
                    x = Gammas,
                    y = all.lambdas,
                    z = log(all.CV.errors)
                ) %>%
                    add_surface() %>%
                    layout(scene = list(
                        xaxis = list(title = "gamma"),
                        yaxis = list(title = "lambda"),
                        zaxis = list(title = "Log CV-Error")
                    ))
            } else if(log.lambda & !log.response){
                plot_ly(
                    x = Gammas,
                    y = log(all.lambdas),
                    z = all.CV.errors
                ) %>%
                    add_surface() %>%
                    layout(scene = list(
                        xaxis = list(title = "gamma"),
                        yaxis = list(title = "log lambda"),
                        zaxis = list(title = "CV-Error")
                    ))
            } else if(log.lambda & log.response){
                plot_ly(
                    x = Gammas,
                    y = log(all.lambdas),
                    z = log(all.CV.errors)
                ) %>%
                    add_surface() %>%
                    layout(scene = list(
                        xaxis = list(title = "gamma"),
                        yaxis = list(title = "log lambda"),
                        zaxis = list(title = "log CV-Error")
                    ))
            }
        })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
