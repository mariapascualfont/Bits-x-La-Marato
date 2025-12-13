df <- read.csv("~/Documents/Clinical cases/CML_df_final.csv")

library(ggplot2)
library(dplyr)

df$resposta_molecular <- factor(
  df$resposta_molecular,
  levels = c("RMP", "RMM", "RM3", "No RM")
)

df$recaiguda <- factor(df$recaiguda, levels = c("no", "sí"))

ggplot(df, aes(x = resposta_molecular)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Distribució de la resposta molecular",
    x = "Resposta molecular",
    y = "Nombre de pacients"
  ) +
  theme_minimal()

ggplot(df, aes(x = recaiguda)) +
  geom_bar(fill = "steelblue") +
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



