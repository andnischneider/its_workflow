ggformat <- theme_classic()+
  theme(axis.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(vjust = 0.5, angle = 45),
        axis.line = element_line(size = 1, linetype = "solid"),
        axis.title = element_text(size = 13, face = "bold"),
        #axis.ticks = element_line(size = 1, linetype = "solid"),
        legend.text = element_text(size=13, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 19, face = "bold"),
        strip.text = element_text(size = 17, face = "bold"))


ggformat_pca <- theme_classic()+
  theme(axis.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(vjust = 1),
        axis.line = element_line(size = 1, linetype = "solid"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.ticks = element_line(size = 1, linetype = "solid"),
        legend.text = element_text(size=13, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 19, face = "bold"),
        strip.text = element_text(size = 17, face = "bold"))

cols_new_treat <- c(rgb(0.30196078431372547, 0.6862745098039216, 0.2901960784313726), 
                    rgb(0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
                    rgb(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))


cols_new_date <- c("#f77189",
              "#97a431",
              "#36ada4",
              "#a48cf4")