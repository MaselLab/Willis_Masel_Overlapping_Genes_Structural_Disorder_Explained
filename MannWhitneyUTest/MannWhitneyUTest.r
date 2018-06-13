library(sm)

# The user should set the working directory by specifying the path to this directory so the example files can be read in.
setwd(dir = "")

# The two example files are read into this script so the relative adaptedness values can be compared
df1 <- scan("./TurnipYellowsVirus_88", sep= "\t")
# The designations file specifies the group to which the relative adaptedness value belongs (either ancestral or novel)
df2 <- scan("./TurnipYellowsVirus_88_Designations", sep="\t")
x_name <- "RelAdapt"
y_name <- "Designation"

# A data frame is then generated using the two files
df = data.frame(df1,df2)
# and a mann-whitney u test is run on the two sets of values. This is where the p-value was found to determine significance
# when comparing ancestral and novel genes. 
wilcox.test(df1~df2, data=df)

# a plot is made for the user so they can visually inspect the codon usage of the two genes
plot(density(df$df1))
group.f <- factor(df2, levels = c(0,1), labels = c("Ancestral", "Novel"))
sm.density.compare(df1,df2,col=c("cyan4","coral"),lwd = 2,xlab = "Relative Adaptedness Values")
colfill <- c(2:(2+length(levels(group.f))))

