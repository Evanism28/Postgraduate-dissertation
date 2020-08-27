# Postgraduate-dissertation
The code manuscript used in this project
library(tidyverse)
library(gridExtra )
library(fPortfolio)
require(graphics)
## Parmeter : Initial Diameter
region = c('rapolano', 'laiatico')
treatment = c('control', 'elevated')
groups = c(paste(region[1], treatment, sep = "_"),paste(region[2], treatment, sep = "_"))
initial_diam = c(1.6, 2.4, 1.7, 3.6) / 1000 # mm to m
as.data.frame(initial_diam) %>% t() -> initial_diam 
colnames(initial_diam) = groups 
initial_diam %>% as.data.frame() -> initial_diam


# File

input_p_rapo = read.csv('P_RA.csv')
input_p_laia = read.csv('P_LA.csv')
dg_rapo = read.csv('RD_RA.csv')
dg_laia = read.csv('RD_LA.csv')
rational_dg = cbind(dg_rapo$d_growth_control, dg_rapo$d_growth_ele, 
                    dg_laia$d_growth_control, dg_laia$d_growth_ele)

rational_dg = rational_dg / 1000 # mm to m

colnames(rational_dg) = groups
rational_dg =  as.data.frame(rational_dg)
rational_dg$year =  dg_laia$year
rational_dg %>% head(1)

input_p_df = cbind(input_p_rapo$p_control, input_p_rapo$p_ele,
                input_p_laia$p_control, input_p_laia$p_ele) %>% as.data.frame()

colnames(input_p_df) = groups
input_p_df$year = dg_laia$year
# input_p_list = cbind(input_p_rapo$p_control, input_p_rapo$p_ele,
#       input_p_laia$p_control, input_p_laia$p_ele) %>% as.data.frame()



## Function
growth_diameter<-function(D, p, rr_star, rs_star, zta, tr){
  P<-p
  As<-(pi*D^2)/4
  a<-118.62
  c<-199.35
  rho<-424
  L<-3.75
  sigma<-10.4
  thof<-1.5
  thor<-tr
  k<-0.72
  y<-0.6
  zeta<- zta
  rr<- P * rr_star
  rs<- P *rs_star
  Hm<-23.31
  H<-Hm*(1-exp(-a*D/Hm))
  z<- H-H^2/a/D
  Ac<-pi*c*D*H/4/a 
  fc<-H/(a*D)
  denominator<-pi/8*rho*(D^2*a*exp(-a*D/Hm)+2*D*H)+L*pi*c/(4*a)*(a*D*(1-H/Hm)+H)*(1/sigma+zeta)
  return(Ac/denominator*(y*(P*(1-exp(-k*L))-rho*(1-fc/2)*H*rs/c-zeta*L*rr)-L*(1/thof/sigma+zeta/thor)))
}

## optimizing
optimizingFun =  function(rrrs = rrrs){
  e = c()
  initial_diam = c(1.6, 2.4, 1.7, 3.6) / 1000
  initial_d = initial_diam[1]#initial_diam$rapolano_control
  p_list = input_p_df$rapolano_control
  dg_acttual = rational_dg$rapolano_control
  d_list = c(initial_d)
  dg_list = c()
  for(i in c(1:length(p_list))){
    d_g = growth_diameter(D = d_list[i], p_list[i], rs_star = rrrs[2], rr_star = rrrs[1], zta = rrrs[3],tr = rrrs[4])
    d_list = c(d_list, d_list[i] + d_g)
    dg_list = c(dg_list, d_g)
  }
  
  a = sum(dg_list < 0) * 1000
  c = ((dg_list - dg_acttual)) %>% sum()
  b = - cor(dg_list, dg_acttual)
  e = c(e, abs(dg_list - dg_acttual) %>% max())
  
  ###########
  
  initial_d = initial_diam[2]#initial_diam$rapolano_control
  p_list = input_p_df$rapolano_elevated
  dg_acttual = rational_dg$rapolano_elevated
  d_list = c(initial_d)
  dg_list = c()
  for(i in c(1:length(p_list))){
    d_g = growth_diameter(D = d_list[i], p_list[i], rs_star = rrrs[2], rr_star = rrrs[1], zta = rrrs[3],tr = rrrs[4])
    d_list = c(d_list, d_list[i] + d_g)
    dg_list = c(dg_list, d_g)
  }
 
  a = a +  sum(dg_list < 0) * 1000
  c = c + ((dg_list - dg_acttual)) %>% sum()
  b = b - cor(dg_list, dg_acttual)
 # e = c(e, abs(dg_list - dg_acttual)) %>% sum()
  e = c(e, abs(dg_list - dg_acttual) %>% max())
  
  #######################
  
  initial_d = initial_diam[3]#initial_diam$rapolano_control
  p_list = input_p_df$laiatico_control
  dg_acttual = rational_dg$laiatico_control
  d_list = c(initial_d)
  dg_list = c()
  for(i in c(1:length(p_list))){
    d_g = growth_diameter(D = d_list[i], p_list[i], rs_star = rrrs[2], rr_star = rrrs[1], zta = rrrs[3],tr = rrrs[4])
    d_list = c(d_list, d_list[i] + d_g)
    dg_list = c(dg_list, d_g)
  }
  
  a = a +  sum(dg_list < 0) * 1000
  c = c + ((dg_list - dg_acttual)) %>% sum()
  b = b - cor(dg_list, dg_acttual)
  #e = c(e, abs(dg_list - dg_acttual)) %>% sum()
  e = c(e, abs(dg_list - dg_acttual) %>% max())
  
  #############################
  initial_d = initial_diam[4]#initial_diam$rapolano_control
  p_list = input_p_df$laiatico_elevated
  dg_acttual = rational_dg$laiatico_elevated
  d_list = c(initial_d)
  dg_list = c()
  for(i in c(1:length(p_list))){
    d_g = growth_diameter(D = d_list[i], p_list[i], rs_star = rrrs[2], rr_star = rrrs[1], zta = rrrs[3],tr = rrrs[4])
    d_list = c(d_list, d_list[i] + d_g)
    dg_list = c(dg_list, d_g)
  }
  
  a = a +  sum(dg_list < 0) * 1000
  c = c + ((dg_list - dg_acttual)) %>% sum()
  b = b - cor(dg_list, dg_acttual)
  e = c(e, abs(dg_list - dg_acttual) %>% max())
  #e = c(e, abs(dg_list - dg_acttual)) %>% sum()
  d = 0.0001 * b +c*10000 /4# +  e *10000/4 
  f = sum(dg_list)
  if(a > 0){
    return(a)
  }else{
    
    return(abs(c) )#b +  c * e *1000 #* -(c + mean(e) +max(e))
  }
}



opt_obj = optim(c(1,0.01, 0.07, 1), optimizingFun,method = 'L-BFGS-B',lower= c(0.5,0.066,0.03,0.45), upper=c(1.5,0.195, 0.105, 1.35) )
#method='L-BFGS-B', lower= c(0.1,0.1), upper=c(0.9,0.9)
#method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#         "Brent"),

rrrs_star = opt_obj$par 

## predict series
growth_seires = function(initial_d, p_list, rs_star, rr_star, zta,tr, dg_acttual = NA){
  d_list = c(initial_d)
  dg_list = c()
  for(i in c(1:length(p_list))){
    d_g = growth_diameter(D = d_list[i], p_list[i], rs_star = rs_star,tr = tr, rr_star = rr_star, zta = zta)
    d_list = c(d_list, d_list[i] + d_g)
    dg_list = c(dg_list, d_g)
  }

  a = sum(dg_list < 0) * 1000
  b = abs(dg_list - dg_acttual) %>% sum()
  # if(a > 0){
  #   return(b)
  # }else{
  #   return(b)
  # }
  return(list(diameter_list =  d_list, diameter_growht_list =  dg_list))
}#3.100275 -1.032573



predict_r_c = growth_seires(initial_d = initial_diam$rapolano_control, 
                             p_list = input_p_df$rapolano_control,
                             rr_star = rrrs_star[1]    , #-1.35352236 ,
                             rs_star = rrrs_star[2] ,#0.05089834 ,
                            zta = rrrs_star[3],
                            tr = rrrs_star[4],
                             dg_acttual = rational_dg$rapolano_control)
rational_dg$rapolano_control_pred = predict_r_c$diameter_growht_list 
rational_dg$rapolano_control_diam = predict_r_c$diameter_list[2:31]

predict_r_e = growth_seires(initial_d = initial_diam$rapolano_elevated, 
                            p_list = input_p_df$rapolano_elevated,
                            rr_star = rrrs_star[1]    , #-1.35352236 ,
                            rs_star = rrrs_star[2] ,#0.05089834 ,
                            zta = rrrs_star[3],
                            tr = rrrs_star[4],
                            dg_acttual = rational_dg$rapolano_elevated)
rational_dg$rapolano_elevated_pred = predict_r_e$diameter_growht_list 
rational_dg$rapolano_elevated_diam = predict_r_e$diameter_list[2:31]


predict_l_c = growth_seires(initial_d = initial_diam$laiatico_control, 
                            p_list = input_p_df$laiatico_control,
                            rr_star = rrrs_star[1]    , #-1.35352236 ,
                            rs_star = rrrs_star[2] ,#0.05089834 ,
                            zta = rrrs_star[3],
                            tr = rrrs_star[4],
                            dg_acttual = rational_dg$laiatico_control)

rational_dg$laiatico_control_pred = predict_l_c$diameter_growht_list 
rational_dg$laiatico_control_diam = predict_l_c$diameter_list[2:31]


predict_l_e = growth_seires(initial_d = initial_diam$laiatico_elevated, 
                            p_list = input_p_df$laiatico_elevated,
                            rr_star = rrrs_star[1]    , #-1.35352236 ,
                            rs_star = rrrs_star[2] ,#0.05089834 ,
                            zta = rrrs_star[3],
                            tr = rrrs_star[4],
                            dg_acttual = rational_dg$laiatico_elevated)


rational_dg$laiatico_elevated_pred = predict_l_e$diameter_growht_list 
rational_dg$laiatico_elevated_diam = predict_l_e$diameter_list[2:31]





# new_groups = paste(groups, "pred", sep =  '_')
# new_groups = c(groups, new_groups) %>% sort()
# rational_dg = rational_dg[,c('year',new_groups)]

library(openxlsx)
write.xlsx(rational_dg, "result_data.xlsx")



# options(repr.plot.width = 14, repr.plot.height = 6)
# df <- rational_dg %>%
#   select(year, rapolano_control_pred, rapolano_control) %>%
#   gather(key = "variable", value = "value", -year)
# p = ggplot(df, aes(x = year, y = value)) + 
#   geom_line(aes(color = variable, linetype = variable)) + 
#   scale_color_manual(values = c("darkred", "steelblue"))
# 
# grid.arrange(p, nrow = 1)

############### Ploting


c('a_sensitive','c_sensitive','Hm_sensitive','rho_sensitive',
  'rr_sensitiv', 'rs_sensitive', 'thof_sensitive', 'thor_sensitive', 'sigma_sensitive', 'L_sensitive', 'zeta_sensitive')



growth_diameter_sensitivey<-function(p= P,D=3, sens){
  P<-p
  As<-(pi*D^2)/4
  a<-118.62
  a = a * sens[1]
  c<-199.35
  c =  c*sens[2]
  rho<-424
  rho =  rho * sens[4]
  L<-3.75
  L = L * sens[10]
  sigma<-10.4
  sigma = sigma * sens[9]
  thof<-1.5
  thof =  thof * sens[7]
  thor<- rrrs_star[4]
  thor = thor * sens[8]
  k<-0.72
  y<-0.6
  zeta<-rrrs_star[3]
  zeta = zeta * sens[11]
  rr<- P * rrrs_star[1]
  rr = rr * sens[5]
  rs<- P * rrrs_star[2]
  rs = rs * sens[6]
  Hm<-23.31
  Hm = Hm * sens[3]
  H<-Hm*(1-exp(-a*D/Hm))
  z<- H-H^2/a/D
  Ac<-pi*c*D*H/4/a 
  fc<-H/(a*D)
  denominator<-pi/8*rho*(D^2*a*exp(-a*D/Hm)+2*D*H)+L*pi*c/(4*a)*(a*D*(1-H/Hm)+H)*(1/sigma+zeta)
  return(Ac/denominator*(y*(P*(1-exp(-k*L))-rho*(1-fc/2)*H*rs/c-zeta*L*rr)-L*(1/thof/sigma+zeta/thor)))
}



seires_sensitivity = function(initial_d, p_list, sens, dg_acttual = NA){
  d_list = c(initial_d)
  dg_list = c()
  for(i in c(1:length(p_list))){
    d_g = growth_diameter_sensitivey(p =p_list[i], D = d_list[i], sens = sens )
    d_list = c(d_list, d_list[i] + d_g)
    dg_list = c(dg_list, d_g)
  }

  return(list(diameter_list =  d_list, diameter_growht_list =  dg_list))
}#3.100275 -1.032573



# FOTO1
sentv_terms = c('a_sensitive','c_sensitive','Hm_sensitive','rho_sensitive',
               'rr_sensitiv', 'rs_sensitive', 'thof_sensitive', 'thor_sensitive', 'sigma_sensitive', 'L_sensitive', 'zeta_sensitive')

sens_level = c(0.5, 1, 1.5)
sentv_matx = matrix(1, nrow = length(sentv_terms) * length(sens_level), ncol = length(sentv_terms))
ini = 1:3
for(i in 1:ncol(sentv_matx)){
  sentv_matx[ini, i] = sens_level
  ini = ini + 3
  
}
for(term in names(initial_diam)){
  graphics.off()
  pdf(paste(term,'_simualtedYRS','.pdf', sep = ''))
  
  par("mar")
  par(mar = c(1, 1, 1, 1), oma = c(4, 4, 0.5, 0.5))
  
  par(mfrow =c(4,3))
  k =  1: 3
  for(t in 1:length(sentv_terms)){
    
    
    windows.options(width=30, height=30)
    tmp_df = matrix(NA, ncol  = 3, nrow = 500) %>% as.data.frame()
    colnames(tmp_df) = paste('scale',sens_level ,sep = "_")
    id = 1
    for(i in k){
      print(k)
      
      tmp_df[, id%%33] = seires_sensitivity(initial_d = initial_diam[[term]],
                                            p_list = input_p_df[[term]] %>% mean() %>% rep(500),
                                            sens = sentv_matx[i,],
                                            dg_acttual = rational_dg[[term]])$diameter_growht_list * 500
      id = id + 1
      
    }
    zoo::as.zoo(tmp_df) %>% ts.plot(col = c('red', 'black', 'blue'), xlab ='ring width(mm)', ylab = 'year')
    legend('bottomright',c((sentv_terms[t] %>% strsplit('_'))[[1]][1]), bty = "n" , cex = 0.8 )
    rm(tmp_df)
    k =  k +3
    
  }
  plot.new()
  legend('center',legend = sens_level, lty = c(1),  col = c('red', 'black', 'blue'),bty = "n" , cex = 1.5 )
  
  mtext("Simulated Year", side = 1, outer = TRUE, cex = 0.9, line = 2.2,  col = "grey20")
  mtext("Ring Width (mm)", side = 2, outer = TRUE, cex = 0.9, line = 2.2,  col = "grey20")
  
  
  
  dev.off()
  
  
}


perc_add = read.csv('Precipitation.csv')
colnames(perc_add) = c('year', 'Laiatioco', 'Rapolano')

### FOTO2
graphics.off()
pdf('ts_dg_two_area.pdf')
par(mfrow = c(2, 1))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 4))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
plot(col = 'red',rational_dg$year,rational_dg$rapolano_elevated_pred *500 ,  xaxt = 'n',type ='o',lty = 1, ylab = c(0.4, 2) )
lines(col ='black',rational_dg$year,rational_dg$rapolano_control_pred * 500, xaxt = 'n',type ='o',lty =2, pch = '*', cex = 1.8)
par(new = T)
plot(perc_add$year, perc_add$Rapolano, col =  'blue',type ='o',lty =4, axes = F, pch = '+', cex = 1.8, xlab=NA, ylab=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Rainfall (mm, total April + May)')


legend('topright', c('elevated CO2', 'control', "rainfall"), lty = c(1,2,4), col = c('red', 'black', 'blue'))
mtext('Rapolano', side = 3, line = -1.5, adj = 0.02, cex = 0.8)

plot(col = 'red',rational_dg$year,rational_dg$laiatico_elevated_pred *500 ,type ='o',lty = 1, ylab = c(0.4, 2))
lines(col = 'black', rational_dg$year,rational_dg$laiatico_control_pred *500, xaxt = 'n',type = 'o', pch ="*",lty =2, cex = 1.8)
par(new = T)
plot(perc_add$year, perc_add$Laiatioco,col ='blue',type ='o',lty =4, axes = F, pch = '+', cex = 1.8, xlab=NA, ylab=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Rainfall (mm, total April + May)')

mtext('Laiatico', side = 3, line = -1.5, adj = 0.02, cex = 0.8)
legend('topright', c('elevated CO2', 'control',"rainfall"), lty = c(1,2,4), col = c('red', 'black', 'blue'))

mtext("Years", side = 1, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
mtext("Ring Width (mm)", side = 2, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
dev.off()



####FOTO21

graphics.off()
pdf('ts_diam_two_area.pdf')
par(mfrow = c(2, 1))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 4))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
plot(rational_dg$year,rational_dg$rapolano_control_diam *1000 ,col = 'red',  xaxt = 'n',type ='o',lty = 1, ylim = c(0, 50) )
lines(rational_dg$year,rational_dg$rapolano_elevated_diam * 1000, col = 'blue',xaxt = 'n',type ='o',lty =2, pch = '*', cex = 1.8)
# par(new = T)
# plot(perc_add$year, perc_add$Rapolano, type ='o',lty =4, axes = F, pch = '+', cex = 1.8, xlab=NA, ylab=NA)
# axis(side = 4)
# mtext(side = 4, line = 3, 'Precipitation')


legend('bottomright', c('elevated CO2', 'control'),col =  c('blue', 'red'),lty = c(1,2))
mtext('Rapolano', side = 3, line = -1.5, adj = 0.02, cex = 0.8)

plot(rational_dg$year,rational_dg$laiatico_control_diam*1000 ,type ='o',col = 'red',lty = 1, ylim = c(0, 50))
lines(rational_dg$year,rational_dg$laiatico_elevated_diam *1000, col = 'blue',xaxt = 'n',type = 'o', pch ="*",lty =2, cex = 1.8)
# par(new = T)
# plot(perc_add$year, perc_add$Laiatioco, type ='o',lty =4, axes = F, pch = '+', cex = 1.8, xlab=NA, ylab=NA)
# axis(side = 4)
# mtext(side = 4, line = 3, 'Precipitation')

mtext('Laiatico', side = 3, line = -1.5, adj = 0.02, cex = 0.8)
legend('bottomright', c('elevated CO2', 'control'), lty = c(1,2), col = c('blue','red'))

mtext("Years", side = 1, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
mtext("Tree Diameter (mm)", side = 2, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
dev.off()

# Foto23

graphics.off()
pdf('ts_dg_two_area_enriched.pdf')
par(mfrow = c(2, 1))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 4))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
plot(col = 'red',rational_dg$year,rational_dg$rapolano_elevated_pred *500 ,  xaxt = 'n',type ='o',lty = 1, ylim = c(0.4, 2) )
lines(col ='blue',rational_dg$year,rational_dg$rapolano_elevated * 500, xaxt = 'n',type ='o',lty =2, pch = '*', cex = 1.8)
#par(new = T)
#plot(perc_add$year, perc_add$Rapolano, col =  'blue',type ='o',lty =4, axes = F, pch = '+', cex = 1.8, xlab=NA, ylab=NA)
#axis(side = 4)
#mtext(side = 4, line = 3, 'Rainfall (mm, total April + May)')


legend('topright', c('simulation', 'real value'), lty = c(1,2), col = c('red',  'blue'))
mtext('Rapolano', side = 3, line = -1.5, adj = 0.02, cex = 0.8)

plot(col = 'red',x = rational_dg$year,y =rational_dg$laiatico_elevated_pred *500 ,  type ='o',lty = 1, ylim = c(0.4, 2) )
lines(col ='blue',rational_dg$year,rational_dg$laiatico_elevated_pred * 500, xaxt = 'n',type ='o',lty =2, pch = '*', cex = 1.8)
#par(new = T)
#plot(perc_add$year, perc_add$Rapolano, col =  'blue',type ='o',lty =4, axes = F, pch = '+', cex = 1.8, xlab=NA, ylab=NA)
#axis(side = 4)
#mtext(side = 4, line = 3, 'Rainfall (mm, total April + May)')


legend('topright', c('simulation', 'real value'), lty = c(1,2), col = c('red',  'blue'))
mtext('Laiatico', side = 3, line = -1.5, adj = 0.02, cex = 0.8)


mtext("Years", side = 1, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
mtext("Ring Width (mm)", side = 2, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
dev.off()


### FOTO3
opt_df = rational_dg
opt_df[2:ncol(opt_df)] = opt_df[2:ncol(opt_df)] * 1000
opt_df[2: 9 ] = opt_df[2: 9 ] /2


# opt_df$log_lai = log(opt_df$laiatico_control_pred)
# opt_df$log_rap = log(opt_df$rapolano_control_pred)
opt_df$log_lai = 1/ opt_df$laiatico_control_pred
opt_df$log_rap = 1 / opt_df$rapolano_control_pred
opt_df$lai_diff = opt_df$laiatico_elevated_pred -  opt_df$laiatico_control_pred
opt_df$rap_diff = opt_df$rapolano_elevated_pred - opt_df$rapolano_control_pred


ols_lai = lm(lai_diff ~ log_lai, data = opt_df)
ols_rapo = lm(rap_diff ~ log_rap, data =  opt_df)

graphics.off()
pdf('regression.pdf')
par(mfrow = c(1, 2))
par( cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 2))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
plot(opt_df$log_lai, opt_df$lai_diff,  pch = 0, cex = 0.8, xlab ="", ylab = "")
mtext('Laiatico', side = 3, line = -1.5, adj = 0.02, cex = 0.8)

abline(ols_lai)

line1 = 'y = a  + b x'
line2 = ": RS, P = pv"
lai_smary = summary(ols_lai)
lai_a = lai_smary$coefficients[1,1] %>% round(3)
lai_b = lai_smary$coefficients[2,1]%>% round(3)
lai_r = lai_smary$r.squared%>% round(3)
lai_p = lai_smary$coefficients[2, 4] %>% round(3)
line1 %>% gsub('a', lai_a,.) %>% gsub('b', lai_b,.) -> l1
line2 %>% gsub('RS', lai_r, .) %>% gsub('pv', lai_p,.) ->l2
mtext(l1, side = 1, outer = FALSE, cex = 0.6, line = -4.5,adj = 0.75,
      col = "grey20")
mtext(expression(paste(R^2,": 0.151", ", P = 0.034"))  , side = 1, outer = FALSE, cex = 0.6, line = -2.5,adj = 0.75,
      col = "grey20")


plot(opt_df$log_rap, opt_df$rap_diff,  pch = 0, 
     cex = 0.8, xlab ="", ylab = "",
     #ylim = c(0, 0.3),
      yaxt = 'n'
     )
axis(side = 4)
mtext('Rapolano', side = 3, line = -1.5, adj = 0.02, cex = 0.8)
ols_lai = ols_rapo
abline(ols_lai)

lai_smary = summary(ols_lai)
lai_a = lai_smary$coefficients[1,1] %>% round(3)
lai_b = lai_smary$coefficients[2,1]%>% round(3)
lai_r = lai_smary$r.squared%>% round(3)
lai_p = lai_smary$coefficients[2, 4] %>% round(3)
line1 %>% gsub('a', lai_a,.) %>% gsub('b', lai_b,.) -> l1
line2 %>% gsub('RS', lai_r, .) %>% gsub('pv', lai_p,.) ->l2
mtext(l1, side = 1, outer = FALSE, cex = 0.6, line = -4.5,adj = 0.75,
      col = "grey20")
mtext(expression(paste(R^2,": 0.047, P = 0.251"))  , side = 1, outer = FALSE, cex = 0.6, line = -2.5,adj = 0.75,
      col = "grey20")

mtext("Ring Width at Control Site(inverse mm)", side = 1, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
mtext("Difference in Ring Width(mm)", side = 2, outer = TRUE, cex = 0.9, line = 2.2,
      col = "grey20")
dev.off()

pdf('ratio.pdf')
opt_df$rap_ratio = opt_df$rapolano_elevated_pred / opt_df$rapolano_control_pred
opt_df$lai_ratio = opt_df$laiatico_elevated_pred / opt_df$laiatico_control_pred
write.xlsx(opt_df, 'result_mm.xlsx')
ts.plot(ts(opt_df[,c("rap_ratio","lai_ratio")], start=c(1964, 1), end=c(1993, 1), frequency=1), 
        col = c('blue', 'red'), xlab = 'Years', ylab = "Ration of Ring Width")
legend('topright', region %>% toupper(), lty = 1, col = c('blue', 'red'), cex = 0.7)
dev.off()
