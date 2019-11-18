# Need these two linraries: Install if missing
library(tidyverse)
library(quantmod)   # For data access

# FRED codes for US GDP growth & CPI level
# Note downloads current data so the results will differ from the blog post
series = c('A191RO1Q156NBEA', 'CPALTT01USQ661S') 
Growth = getSymbols(series[1], src='FRED', auto.assign=FALSE) 
CPI    = getSymbols(series[2], src='FRED', auto.assign=FALSE) 

Data = inner_join(tibble(Date=time(Growth), Growth=coredata(Growth)), 
                  tibble(Date=time(CPI), CPI=coredata(CPI)), by=c("Date")) %>% 
  mutate(Inflation=100*(CPI/lag(CPI,4)-1)) %>%
  select(Date, Growth, Inflation) %>% 
  drop_na() # Drop missing obs to balance dataset

centre_colour = c("seagreen","tomato") # Colours for time series/centre of fancharts
tail_colour   = "gray95"               # Colour for the tails, used later but defined here 
pivot_longer(Data, cols=-Date, names_to="Variables", values_to="Values") %>% 
  ggplot() + 
  geom_line(aes(x=Date, y=Values, group=Variables, colour=Variables), size=1.1, show.legend=TRUE) +
  scale_colour_manual(values=centre_colour) +
  theme_minimal() + 
  theme(legend.title = element_blank()) +
  labs(title="US GDP growth and CPI inflation", x="", y="",
       caption=paste0("Source: FRED series ", paste(series, collapse=", ")))

# Model

m     = 4  # maximum lag in VAR
Datal = Data %>%
  pivot_longer(cols=-Date, names_to="Names", values_to="Values") %>%
  mutate(lag_value=list(0:m)) %>%
  unnest(cols=lag_value) %>%
  group_by(Names, lag_value) %>%
  mutate(Values=lag(Values, unique(lag_value))) %>%
  ungroup() %>%
  mutate(Names = if_else(lag_value==0, Names,                                      # No suffix at lag 0
                         paste0(Names, "_", str_pad(lag_value, 2, pad="0")))) %>%  # All other lags
  select(-lag_value) %>%      # Drop the redundant lag index
  pivot_wider(names_from=Names, values_from=Values) %>%
  slice(-c(1:m)) %>%          # Remove missing lagged initial values
  mutate(constant = 1)        # Add column of ones at end

s = paste(paste0(str_pad(1:m, 2, pad="0"), "$"), collapse="|")
X = data.matrix(select(Datal,  matches(paste0(s,"|constant"))))
Y = data.matrix(select(Datal, -matches(paste0(s,"|constant|Date"))))

(bhat = solve(crossprod(X), crossprod(X,Y)))

# Forecast

nv    = ncol(Y) # Number of variables
nf    = 12      # Periods to forecast
nb    = 16      # Periods of back data to plot, used later

v     = crossprod(Y - X %*% bhat)/(nrow(Y)-m*nv-1)            # Calculate error variance
bhat2 = bhat[rep(seq(1,m*nv,m),m) + rep(seq(0,m-1), each=nv),] # Reorder for simulation
A     = rbind(t(bhat2), diag(1,nv*(m-1), nv*m))                # First order form - A 
B     = diag(1,nv*m,nv)                                        # First order form - B
cnst  = c(t(tail(bhat,1)), rep(0,nv*(m-1)))                    # First order constants

# Simulation loop
Yf     = matrix(0,nv*m,nf+1)                # Stores forecasts
Yf[,1] = c(t(tail(Y,m)[m:1,]))              # Lagged data
Pf     = matrix(0,nv,nf+1)                  # Stores variances
P      = matrix(0,nv*m,nv*m)                # First period state covariance

for (k in 1:nf) { 
  Yf[,k+1] = cnst + A %*% Yf[,k]
  P        = A %*% P %*% t(A) + B %*% v %*% t(B)
  Pf[,k+1] = diag(P)[1:nv]
  }

# Fancharts
qu     = c(.05,.2,.35,.65,.8,.95)  # Chosen quantiles ensures 30% of the distribution each colour
nq     = length(qu)
fdates = seq.Date(tail(Data$Date,1), by="quarter", length.out=nf+1) # Forecast dates
forecast_data = tibble(Date     = rep(fdates, 2), 
                       Variable = rep(colnames(Data)[-1], each=(nf+1)), 
                       Forecast = c(t(Yf[1:nv,])),
                       Variance = c(t(sqrt(Pf)))) %>% 
  bind_cols(map(qu, qnorm, .$Forecast, .$Variance)) %>%         # Calculate quantiles
  select(-c("Forecast", "Variance")) %>% 
  {bind_rows(select(., -(nq+2)),                                # Drop last quantile 
             select(., -3) %>%                                  # Drop first quantile
               arrange(Variable, desc(Date)) %>%                # Reverse order
               rename_at(-(1:2), ~paste0("V",1:(nq-1))) )} %>%  # Shift names of reversed ones 
  pivot_longer(cols=-c(Date, Variable), names_to="Area", values_to="Coordinates") %>% 
  unite(VarArea, Variable, Area, remove=FALSE) %>%              # Create variable to index polygons
  bind_rows(pivot_longer(tail(Data,nb), cols = -Date, names_to="Variable", values_to="Backdata"), .)

# Band colours 'ramp' from the centre to the tail colour 
band_colours = colorRampPalette(c(rbind(tail_colour, centre_colour), tail_colour), 
                                space="Lab")(nv*nq+1)[-seq(1, nv*nq+1, nq)]

# Final plot
ggplot(forecast_data) + 
  geom_rect(aes(xmin=Date[nv*nb], xmax=max(Date), ymin=-Inf, ymax=Inf), fill=tail_colour, alpha=.2) +  
  geom_polygon(aes(x=Date, y=Coordinates, group=VarArea, fill=VarArea)) +
  scale_fill_manual(values=band_colours) +
  geom_line(aes(x=Date, y=Backdata, group=Variable, colour=Variable)) +
  scale_colour_manual(values=centre_colour) +
  scale_x_date(expand=c(0,0)) +
  theme_minimal() +
  theme(legend.position="none") +
  facet_wrap(~ Variable, ncol=1) +
  labs(title="Forecasts of US GDP growth and CPI inflation", 
       subtitle=paste("Quarterly data, annual rates of change, VAR with", m, "lags"), 
       caption=paste("Source: FRED series", paste(series, collapse=", ")), x="", y="")

