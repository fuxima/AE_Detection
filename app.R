library(tidyverse)
library(survival)
library(survminer)
library(openEBGM)
library(tidytext)
library(lubridate)
library(shiny)
library(DT)

# Define UI
ui <- fluidPage(
  plotOutput("ae_frequency"),
  plotOutput("survival_plot"),
  DTOutput("prr_table")
)

# Define server logic
server <- function(input, output) {
  # Create reactive data and analysis
  analysis_results <- reactive({
    # 1. Create sample data
    set.seed(123)
    ae_data <- tibble(
      patient_id = rep(1:100, each = 2),
      term = sample(c("Headache", "Nausea", "Rash", "Fatigue", "Dizziness", "Fever"), 
                    200, replace = TRUE, 
                    prob = c(0.3, 0.2, 0.15, 0.15, 0.1, 0.1)),
      severity = sample(c("Mild", "Moderate", "Severe"), 200, replace = TRUE),
      onset_date = as.Date("2023-01-01") + days(sample(0:90, 200, replace = TRUE)),
      drug = sample(c("DrugA", "DrugB"), 200, replace = TRUE)
    )
    
    clinical_notes <- tibble(
      note_id = 1:50,
      text = paste(
        "Patient reported", 
        sample(c("headache", "nausea", "rash", "fatigue", "dizziness", "fever"), 50, replace = TRUE),
        sample(c("after medication", "during treatment", "at follow-up"), 50, replace = TRUE)
      )
    )
    
    # Return both datasets and results in a list
    list(
      ae_data = ae_data,
      clinical_notes = clinical_notes,
      results = {
        # 2. Analysis functions
        compute_prr <- function(data, drug_name, ae_term) {
          a <- sum(data$drug == drug_name & data$term == ae_term)
          b <- sum(data$drug != drug_name & data$term == ae_term)
          c <- sum(data$drug == drug_name & data$term != ae_term)
          d <- sum(data$drug != drug_name & data$term != ae_term)
          
          (a/(a+c)) / (b/(b+d))
        }
        
        ae_pipeline <- function(data, notes) {
          # Frequency analysis
          frequent_events <- data %>% 
            count(term, severity, sort = TRUE)
          
          # Time-to-event analysis
          time_fit <- survfit(Surv(as.numeric(onset_date - min(onset_date)), rep(1, nrow(data))) ~ drug, data)
          
          # PRR calculation
          prr_table <- cross_df(list(drug = unique(data$drug), 
                                     term = unique(data$term))) %>%
            mutate(prr = map2_dbl(drug, term, ~compute_prr(data, .x, .y)))
          
          # NLP extraction
          nlp_results <- notes %>%
            unnest_tokens(word, text) %>%
            semi_join(data, by = c("word" = "term"))
          
          list(
            frequent_events = frequent_events,
            time_fit = time_fit,
            prr_table = prr_table,
            nlp_results = nlp_results
          )
        }
        
        # 3. Run analysis pipeline
        ae_pipeline(ae_data, clinical_notes)
      }
    )
  })
  
  # Outputs
  output$ae_frequency <- renderPlot({
    all_data <- analysis_results()
    ggplot(all_data$results$frequent_events, aes(x = reorder(term, n), y = n, fill = severity)) +
      geom_col() + coord_flip() +
      labs(title = "Adverse Event Frequency", x = "AE Term", y = "Count")
  })
  
  output$survival_plot <- renderPlot({
    all_data <- analysis_results()
    ggsurvplot(all_data$results$time_fit, data = all_data$ae_data, risk.table = TRUE)$plot
  })
  
  output$prr_table <- renderDT({
    all_data <- analysis_results()
    DT::datatable(all_data$results$prr_table %>% arrange(desc(prr))) %>%
      DT::formatRound("prr", 2) %>%
      DT::formatStyle("prr", 
                      background = DT::styleColorBar(range(all_data$results$prr_table$prr), 'lightblue'),
                      backgroundSize = '98% 88%')
  })
}

# Run the application
shinyApp(ui = ui, server = server)