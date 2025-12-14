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

df$Tipus.transcript <- factor(
  df$Tipus.transcript,
  levels = c('e14a2', 'e13a2', 'e1a2', 'Unknown')
)

ggplot(df, aes(x = Tipus.transcript)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  labs(
    title = "Distribució de Tipus de Transcript",
    x = "Tipus de Transcript",
    y = "Nombre de pacients"
  ) +
  theme_minimal() 

ggplot(
  df %>% filter(!is.na(Tipus.transcript), !is.na(Grup.Edat)),
  aes(x = Tipus.transcript, fill = Grup.Edat)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de Tipus de Transcript segons edat",
    x = "Tipus de Transcript",
    y = "Proporció",
    fill = "Edat"
  ) +
  theme_minimal()

ggplot(
  df %>% filter(!is.na(Recaiguda), !is.na(Tipus.transcript)),
  aes(x = Recaiguda, fill = Tipus.transcript)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de Recaiguda segons Tipus de Transcript",
    x = "Recaiguda",
    y = "Proporció",
    fill = "Tipus de Transcript"
  ) +
  theme_minimal()

ggplot(
  df %>% filter(!is.na(Resposta.Molecular), !is.na(Tipus.transcript)),
  aes(x = Resposta.Molecular, fill = Tipus.transcript)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proporció de Resposta Molecular segons Tipus de Transcript",
    x = "Resposta Molecular",
    y = "Proporció",
    fill = "Tipus de Transcript"
  ) +
  theme_minimal()
