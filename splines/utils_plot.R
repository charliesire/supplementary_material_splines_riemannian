
arrange_plots = function(plotlist, ncol, legend, title){ #function to plot multiplot ggplot 
  res_plot = plot_grid(plotlist = plotlist, ncol=ncol)
  res_plot = plot_grid(res_plot,legend, ncol=2, rel_widths = c(10,1))
  final_plot <- ggdraw() + 
    draw_plot(res_plot) + 
    draw_label(title, size = 12, fontface = 'bold', x = 0.5, y = 1.1, hjust = 0.5) + theme(plot.margin = margin(t = 60)) 
}
# Extract legend


plot_multiple_df = function(list_df, coords, value, name = "", list_titles = "", big_title = "", list_opt = NULL, ncol = 4, all_col = FALSE,add_func = NULL){
  if(length(value)>1){
    for (xx in 1:length(value)){colnames(list_df[[xx]])[which(colnames(list_df[[xx]]) == value[xx])] = "aa"}
    value = "aa"
  }
  values_for_col =  unlist(lapply(list_df, function(x){x[,value]}))
  values_for_col=values_for_col[!is.na(values_for_col)]
  list_plot = list()
  if(all_col){colors = c("purple", "blue", "green", "white" ,"yellow","orange", "red")}
  else{colors = magma(256)}
  for(idx in seq_len(length(list_df))){
    width = max(diff(as.numeric(as.matrix(list_df[[idx]][coords[1]]))))
    list_plot[[idx]] = ggplot() + geom_tile(data = list_df[[idx]],  aes_string(x = coords[1],y = coords[2], fill = value, width = width))+
      scale_fill_gradientn(colors = colors,limits = c(min(values_for_col), max(values_for_col)), name = name) + theme_minimal() + labs(title = list_titles[idx]) + theme(plot.title = element_text(size = 20, face = "bold")) 
    if(!is.null(list_opt)){
      colnames(list_opt[[idx]]) = c("z","theta")
      list_plot[[idx]] = list_plot[[idx]] + geom_point(data = list_opt[[idx]], aes(x = z, y = theta), color = "black", size = 1, alpha = 1)
    }
    if(!is.null(add_func)){list_plot[[idx]] = list_plot[[idx]] + add_func()}
    if(idx == 1){legend = get_legend(list_plot[[idx]]+ theme(legend.key.size = unit(0.6, 'cm'), legend.title=element_blank(), legend.text = element_text(size=12.5)))}
    list_plot[[idx]] = list_plot[[idx]] + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  res = arrange_plots(plotlist = list_plot, ncol = ncol, legend = legend, title = big_title)
  return(res)
}

plot_multiple_columns = function(df, coords, list_value, name = "", list_titles = "", big_title = "", list_opt = NULL, ncol = 4, all_col = FALSE, add_func = NULL, remove_na = FALSE){
  df = df[complete.cases(df), ]
  values_for_col =  unlist(sapply(list_value, function(x){df[,x]}))
  values_for_col=values_for_col[!is.na(values_for_col)]
  list_plot = list()
  if(all_col){colors = colors = c("purple", "blue", "green", "white" ,"yellow","orange", "red")}
  else{colors = magma(256)}
  for(idx in seq_len(length(list_value))){
    list_plot[[idx]] = ggplot() + geom_tile(data = df,  aes_string(x = coords[1],y = coords[2], fill = list_value[idx]))+
      scale_fill_gradientn(colors = colors,limits = c(min(values_for_col), max(values_for_col)), name = name) + theme_minimal() + labs(title = list_titles[idx]) + theme(plot.title = element_text(size = 10, face = "bold")) 
    if(!is.null(list_opt)){
      colnames(list_opt[[idx]]) = c("z","theta")
      list_plot[[idx]] = list_plot[[idx]] + geom_point(data = list_opt[[idx]], aes(x = z, y = theta), color = "black", size = 1, alpha = 1)
    }
    if(!is.null(add_func)){list_plot[[idx]] = list_plot[[idx]] + add_func()}
    if(idx == 1){legend = get_legend(list_plot[[idx]]+ theme(legend.key.size = unit(0.6, 'cm'), legend.title=element_blank(), legend.text = element_text(size=12.5)))}
    list_plot[[idx]] = list_plot[[idx]] + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  res = arrange_plots(plotlist = list_plot, ncol = ncol, legend = legend, title = big_title)
  return(res)
}