library(tuneR)
library(shiny)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(janitor)
library(knitr)
normalify <- function(seq){
  return(qnorm((2*rank(seq)-1)/(2*length(seq))))
}
permcov <- function(x,y,nsims = 10000){
  n = length(x)
  if (n < 2){return(c(0,0))}
  sample_cov = cov(x,y)
  more_cov = 0
  for (i in 1:nsims){
    perm = sample(c(1:n),n)
    cur_t = cov(x,y[perm])
    if (abs(cur_t) > abs(sample_cov)){
      more_cov = more_cov + 1
    } 
  }
  return(c(sample_cov, more_cov/nsims))
}
decimals <- function(x,n){
  return(format(round(x, n), nsmall=n))
}
covid_raw <- read.csv("owid-covid-data.csv") %>% janitor::clean_names()
covid_cleaner <- covid_raw
covid_cleaner <- covid_cleaner[]
covid_cleaner[covid_cleaner=="NaN"] <- NA
covid_cleaner[covid_cleaner==-Inf] <- NA
#A vector containing the names of locations which are greater than a single nation. 
super_locations = c("Africa", "Asia", "European Union", "High income", "International", "Low income", "Lower middle income", "North America", "Oceania", "South America", "Upper middle income", "World")
covid_cleaner <- covid_cleaner[!(covid_cleaner$location %in% super_locations),]

covid = covid_cleaner
covid <- covid %>% group_by(location)
covid = covid %>% summarise(
  continent = max(continent),
  total_vaccinations_per_hundred = max(total_vaccinations_per_hundred, na.rm = T),
  people_vaccinated_per_hundred = max(people_vaccinated_per_hundred, na.rm = T),
  people_fully_vaccinated_per_hundred = max(people_fully_vaccinated_per_hundred, na.rm = T),
  total_boosters_per_hundred = max(total_boosters_per_hundred, na.rm = T),
  gdp_per_capita = mean(gdp_per_capita, na.rm=T),
  log_gdp_per_capita = log(mean(gdp_per_capita, na.rm=T),base=10),
  population = mean(population, na.rm = T)
)
covid[covid=="NaN"] <- NA
covid[covid==-Inf] <- NA
covid <- covid[!is.na(covid$gdp_per_capita),]

#"Asia"          "Europe"        "Africa"        "North America" "South America" "Oceania"
ui <- fluidPage(
  titlePanel("Comparing vaccination rates with gdp per capita"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "vax",
                  label = "Vaccination metric:",
                  choices = list("Vaccinations per capita" = "total_vaccinations_per_hundred",
                                 "People vaccinated per capita" = "people_vaccinated_per_hundred",
                                 "People fully vaccinated per capita"="people_fully_vaccinated_per_hundred" ,
                                 "Boosters per capita"="total_boosters_per_hundred")
      ),
      checkboxInput(inputId = "contas",
                    label = "Include locations in Asia",
                    value = T,
                    width = NULL),
      checkboxInput(inputId = "conteu",
                    label = "Include locations in Europe",
                    value = T,
                    width = NULL),
      checkboxInput(inputId = "contaf",
                    label = "Include locations in Africa",
                    value = T,
                    width = NULL),
      checkboxInput(inputId = "contna",
                    label = "Include locations in North America",
                    value = T,
                    width = NULL),
      checkboxInput(inputId = "contsa",
                    label = "Include locations in South America",
                    value = T,
                    width = NULL),
      checkboxInput(inputId = "contoc",
                    label = "Include locations in Oceania",
                    value = T,
                    width = NULL),
      checkboxInput(inputId = "combine",
                    label = "Use range union (instead of intersection)",
                    value = FALSE,
                    width = NULL),
      checkboxInput(inputId = "ranked",
                    label = "Ranked values",
                    value = FALSE,
                    width = NULL),
      sliderInput(inputId = "poprange",label = "Range of locations included by population percentile:",
                  min = 0,
                  max = 100,
                  value = c(0,100),
                  sep = ""),
      sliderInput(inputId = "gdprange",
                  label = "Range of locations included by gdp per capita percentile:",
                  min = 0,
                  max = 100,
                  value = c(0,100),
                  sep = ""),
      submitButton(text = "Create plot")
    ),
    mainPanel(
      plotOutput(outputId = "nameplot"),
      verbatimTextOutput("verb")
    )
  )
)
server <- function(input,output) {
  
  
  
  output$nameplot <- renderPlot({
    included_continents = c("Asia","Europe","Africa","North America","South America","Oceania")[c(
      input$contas,input$conteu,input$contaf,input$contna,input$contsa,input$contoc
    )]
    
    
    
    
    contrel = covid$continent %in% included_continents
    names = covid$location[contrel]
    pop = covid$population[contrel]
    gdp = covid$gdp_per_capita[contrel]
    n = length(gdp)
    gdprel = order(gdp)[(floor(n-1)*input$gdprange[1]/100+1):ceiling((n-1)*input$gdprange[2]/100+1)]
    poprel = order(pop)[(floor(n-1)*input$poprange[1]/100+1):ceiling((n-1)*input$poprange[2]/100+1)]
    rel = intersect(gdprel,poprel)
    if (input$combine){
      rel = union(gdprel, poprel)
    }
    names = names[rel]
    relevant_rows = covid$location %in% names
    y = pull(covid, input$vax)[relevant_rows]/100
    x = covid$log_gdp_per_capita[relevant_rows]
    continents = covid$continent[relevant_rows]
    pop = covid$population[relevant_rows]
    if (input$ranked){
      y = rank(y)
      x = rank(x)
    }
    data = data.frame(x=x,y=y,labels = names,continent = continents)
    data <- na.omit(data)
    x = data$x
    y = data$y
    yaxislabel = list("total_vaccinations_per_hundred" = "Total vaccinations per capita",
                      "people_vaccinated_per_hundred" = "People vaccinated per capita",
                      "people_fully_vaccinated_per_hundred" = "People fully vaccinated per capita",
                      "total_boosters_per_hundred" = "Total boosters per capita")[input$vax]
    if (input$ranked){
      yaxislabel = list("total_vaccinations_per_hundred" = "Ranked total vaccinations per capita",
                        "people_vaccinated_per_hundred" = "Ranked people vaccinated per capita",
                        "people_fully_vaccinated_per_hundred" = "Ranked people fully vaccinated per capita",
                        "total_boosters_per_hundred" = "Ranked total boosters per capita")[input$vax]
    }
    
    
    xaxislabel = "Log gdp per capita"
    if (input$ranked){
      xaxislabel = "Ranked gdp per capita"
    }
    
    results = permcov(x,y) 
    m = results[1]/var(x)
    b = mean(data$y) - m * mean(data$x)
    p = results[2]
    
    msg1 = paste("This data has a correlation of '",decimals(cor(x,y), 4),"' and a p-value of '",decimals(p, 3),"'.\n    (p-value from two-sided permutation test with H0 of independance)",sep='')
    msg2 = paste("Additionally, using linear regression produces the model: \n   '",yaxislabel,"' = '",xaxislabel,"' * ",decimals(m,3)," + ",decimals(b,3),sep='')
    output_message = paste(msg1,'\n',msg2)
    output$verb <- renderText({ output_message})
    
    return(data %>%
             ggplot(aes(x = x,
                        y = y)) +
             geom_point(aes(colour = continent)) + 
             theme_minimal() + 
             xlab(xaxislabel) +
             ylab(yaxislabel) +
             geom_smooth(method = "lm", se = FALSE)+
             theme(axis.text=element_text(size=12),
                   axis.title=element_text(size=14,face="bold"))
    )
  })
}
shinyApp(ui = ui, server = server)
