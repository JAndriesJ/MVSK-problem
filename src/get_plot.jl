module get_plot
## Plotting the data
plot(df.time, log.([df.DAC df.IBM df.SAVA df.TSCO]),
xlabel="time",
ylabel="log. close price",
label=["DAC" "IBM" "SAVA" "TSCO"],
xrotation=70)

## Computing Empirical moments

time(df.time[1])

mean([df.DAC df.IBM df.SAVA df.SHOP df.TSCO],dims=1)
log.([df.DAC df.IBM df.SAVA df.SHOP df.TSCO])



end