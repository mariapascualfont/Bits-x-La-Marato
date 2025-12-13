df <- read.csv("~/Documents/Clinical cases/CML_df_final.csv")

library(ggplot2)
library(dplyr)

df$resposta_molecular <- factor(
  df$resposta_molecular,
  levels = c("RMP", "RMM", "RM3", "No RM")
)

df$recaiguda <- factor(df$recaiguda, levels = c("no", "sí"))
ggplot(df, aes(x = resposta_molecular)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de la resposta molecular",
    x = "Resposta molecular",
    y = "Nombre de pacients"
  ) +
  theme_minimal()

ggplot(df, aes(x = recaiguda)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de la recaiguda",
    x = "Recaiguda",
    y = "Nombre de pacients"
  ) +
  theme_minimal()


ggplot(df, aes(x = resposta_molecular, fill = recaiguda)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de recaiguda segons resposta molecular",
    x = "Resposta molecular",
    y = "Proporció",
    fill = "Recaiguda"
  ) +
  theme_minimal()

df$SEXO <- factor(
  df$SEXO,
  levels = c("0", "1"), 
  labels = c("Dona", "Home")
)

ggplot(df, aes(x = SEXO)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució del sexe dels pacients",
    x = "Sexe",
    y = "Nombre de pacients"
  ) +
  theme_minimal()

ggplot(df, aes(x = resposta_molecular, fill = SEXO)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de resposta molecular segons sexe",
    x = "Resposta molecular",
    y = "Proporció",
    fill = "Sexe"
  ) +
  theme_minimal()

ggplot(df, aes(x = recaiguda, fill = SEXO)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de recaiguda segons sexe",
    x = "Recaiguda",
    y = "Proporció",
    fill = "Sexe"
  ) +
  theme_minimal()

library(dplyr)
library(lubridate)
library(stringr)
convertir_data <- function(x) {
  
  x <- tolower(str_trim(x))
  
  # Separar dia, mes, any
  parts <- str_split_fixed(x, "-", 3)
  
  dia  <- as.integer(parts[,1])
  mes_txt <- parts[,2]
  any  <- as.integer(parts[,3])
  
  # Convertir any 2 dígits → 4 dígits
  any <- ifelse(any > 25, 1900 + any, 2000 + any)
  
  # Diccionari mesos
  mesos <- c(
    "ene" = 1, "feb" = 2, "mar" = 3, "abr" = 4,
    "may" = 5, "jun" = 6, "jul" = 7, "ago" = 8,
    "sep" = 9, "oct" = 10, "nov" = 11, "dic" = 12
  )
  
  mes <- mesos[mes_txt]
  
  # Crear Date manualment
  make_date(year = any, month = mes, day = dia)
}


df$data_naixement <- convertir_data(df$Data.naixement)



df <- df %>%
  mutate(
    edat = floor(time_length(interval(data_naixement, Sys.Date()), "years"))
  )

df <- df %>%
  mutate(
    edat_grup = cut(
      edat,
      breaks = c(0, 40, 60, 80, Inf),
      labels = c("<40", "40–59", "60–79", "≥80"),
      right = FALSE
    )
  )

ggplot(df, aes(x = edat_grup)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de l'edat",
    x = "Grup d'edat",
    y = "Nombre de pacients"
  ) +
  theme_minimal() 

ggplot(
  df %>% filter(!is.na(edat_grup), !is.na(resposta_molecular)),
  aes(x = resposta_molecular, fill = edat_grup)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de resposta molecular segons edat",
    x = "Resposta molecular",
    y = "Proporció",
    fill = "Edat"
  ) +
  theme_minimal()

ggplot(
  df %>% filter(!is.na(edat_grup), !is.na(recaiguda)),
  aes(x = recaiguda, fill = edat_grup)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de recaiguda segons edat",
    x = "Recaiguda",
    y = "Proporció",
    fill = "Edat"
  ) +
  theme_minimal()
