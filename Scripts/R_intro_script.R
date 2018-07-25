# Introduction to R script!

#### INTRODUCTION ####
# make a character vector of disney characters
characters <- c("Cinderella", "Belle", "Dory")

# check the class of 'characters'
class(characters)

# make a numeric vector of dates
dates <- c(1950, 1987, 2003)

# check the class of 'dates'
class(dates)

# make a logical vector of princess
princess <- c(TRUE, TRUE, FALSE)
    
# check the class of 'princess'
class(princess)
    
# combine 'characters', 'dates', and 'princess' to make a new dataframe
disney <- data.frame(characters, dates, princess)

# view your new dataframe
View(disney)

# look at the structure of the dataframe
str(disney)

# change the characters column to be a character not factor
disney$characters <- as.character(disney$characters)

# how many rows does the dataframe have?
nrow(disney)

# how many columns does the dataframe have?
ncol(disney)

180# looking at your dataframe by row
disney[1,]

# looking at your dataframe by column
disney[,1]

# looking at one 'cell' of your dataframe
disney[1,1]
    
#### LOGIC ####

# which elements of characters (if any) are equal to "Cinderella"?
characters == "Cinderella"

# what elements of characters are not equal to "Cinderella"?
characters != "Cinderella"

# is "Cinderella" in the characters vector?
"Cinderella" %in% characters

# Are there movies older than 1995?
dates < 1995

# IF/ELSE statements - print "Hello, Cinderella" if character is Cinderella, otherwise print "You're not Cinderella"
ifelse(characters == "Cinderella", "Hello, Cinderella", "You're not Cinderella")

if(characters == "Cinderella") {
    print("Hello, Cinderella")
}

# ON YOUR OWN (OYO) - use if/else to print "Before me!" or "After me!" for each movie based on if it came out before or after you were born
ifelse(dates < 1993, "Before me!", "After me!")

# For loops - introduction
for(i in 1:10) {
    print(i)
}

# For loops - print each character
for(row in 1:nrow(disney)) {
    print(disney$character[row])
}

# For loops - more advanced
# print each character IF they are a princess ELSE print "Not a princess."
for(row in 1:nrow(disney)) {
    ifelse(disney$princess[row] == TRUE, print(disney$characters[row]), print("Not a princess."))
}

#### PACKAGES ####

# install dplyr, tidyr, and ggplot2
install.packages(c("dplyr", "tidyr", "ggplot2"))

