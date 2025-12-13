df <- read.csv("~/Documents/Clinical cases/CML_df_final.csv")

library(ggplot2)
library(dplyr)

df$Resposta.Molecular <- factor(
  df$Resposta.Molecular,
  levels = c("RMP", "RMM", "RM3", "No RM")
)

df$Recaiguda <- factor(df$Recaiguda, levels = c("no", "sí"))
ggplot(df, aes(x = Resposta.Molecular)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de la resposta molecular",
    x = "Resposta molecular",
    y = "Nombre de pacients"
  ) +
  theme_minimal()

ggplot(df, aes(x = Recaiguda)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de la Recaiguda",
    x = "Recaiguda",
    y = "Nombre de pacients"
  ) +
  theme_minimal()


ggplot(df, aes(x = Resposta.Molecular, fill = Recaiguda)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de Recaiguda segons resposta molecular",
    x = "Resposta molecular",
    y = "Proporció",
    fill = "Recaiguda"
  ) +
  theme_minimal()

df$Sexe <- factor(
  df$Sexe,
  levels = c("Dona", "Home"), 
)

ggplot(df, aes(x = Sexe)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució del sexe dels pacients",
    x = "Sexe",
    y = "Nombre de pacients"
  ) +
  theme_minimal()

ggplot(df, aes(x = Resposta.Molecular, fill = Sexe)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de resposta molecular segons sexe",
    x = "Resposta molecular",
    y = "Proporció",
    fill = "Sexe"
  ) +
  theme_minimal()

ggplot(df, aes(x = Recaiguda, fill = Sexe)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de Recaiguda segons sexe",
    x = "Recaiguda",
    y = "Proporció",
    fill = "Sexe"
  ) +
  theme_minimal()

df$Grup.Edat <- factor(
  df$Grup.Edat,
  levels = c("<40", "40–59", "60-79", "≥80"), 
)

ggplot(df, aes(x = Grup.Edat)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de l'edat",
    x = "Grup d'edat",
    y = "Nombre de pacients"
  ) +
  theme_minimal() 

ggplot(
  df %>% filter(!is.na(Grup.Edat), !is.na(Resposta.Molecular)),
  aes(x = Resposta.Molecular, fill = Grup.Edat)
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
  df %>% filter(!is.na(Grup.Edat), !is.na(Recaiguda)),
  aes(x = Recaiguda, fill = Grup.Edat)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de Recaiguda segons edat",
    x = "Recaiguda",
    y = "Proporció",
    fill = "Edat"
  ) +
  theme_minimal()
