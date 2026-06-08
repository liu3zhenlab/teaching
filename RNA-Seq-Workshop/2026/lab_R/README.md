## **1. Chi-square test**

```r
d <- c(12, 36, 24, 70)
dm <- matrix(d, nrow=2, byrow=T)
dm
chisq.test(dm)
```

### **What it does**

Create a vector:

```r
d
# 12 36 24 70
```

Convert it into a 2 × 2 matrix:

```r
     [,1] [,2]
[1,]   12   36
[2,]   24   70
```

Perform a chi-square test:

```r
chisq.test(dm)
```

Used to test whether the row and column categories are independent.

------

## **2. Graphics settings (commented out)**

```r
#options(bitmapType = 'cairo')
#Sys.setenv("DISPLAY"=":0.0")
```

These lines are disabled by `#`.

They are sometimes needed on Linux servers to enable plotting.

------

## **3–4. Load custom functions**

```r
ts="https://.../xcard.R"
source(ts)
trend="https://.../trend.R"
source(trend)
```

`source()` downloads and executes an R script.

Equivalent to copying the code into your current session.

------

## **5. Create and enter a working directory**

```r
dir.create("~/rnaseq")
setwd("~/rnaseq")
```

### **Purpose**

Create folder:

```text
~/rnaseq
```

and make it the current working directory.



Check current directory:

```r
getwd()
```

------

## **6. Basic arithmetic**

```r
2 + 4
68 * 0.15
```

Results:

```r
6
10.2
```

------

## **7. Variables**

```r
y <- 2
y = 2
y
```

Both assignment operators work.

```r
info <- "hello world"
cat(info)
```

`cat()` prints text without quotes.



Output:

```text
hello world
```

------

## **8. Assignment example**

```r
y <- 2 + 4
```

Stores 6 in variable `y`.

------

## **9. Case sensitivity**

```r
y <- 2
Y <- 3

y
Y
```

Output:

```r
2
3
```

R treats uppercase and lowercase as different variables.

------

# **Vectors**

## **10. Numeric vectors**

```r
x <- c(10.4, 5.6, 3.1, 6.4, 21.7)
```

Creates a numeric vector.

### **Sum**

```r
sum(x)
```

Result:

```r
47.2
```

### **Multiply every element**

```r
2*x
```

Result:

```r
20.8 11.2 6.2 12.8 43.4
```

### **Extract second element**

```r
x[2]
```

Result:

```r
5.6
```

------

## **11. Logical vectors**

```r
lv <- c(TRUE, FALSE, TRUE, TRUE)
```

### **Negation**

```r
!lv
```

Result:

```r
FALSE TRUE FALSE FALSE
```

### **Compare**

```r
lv == FALSE
```

Result:

```r
FALSE TRUE FALSE FALSE
```

------

## **12. Character vectors**

```r
cv <- c("a", "b", "c")
```

### **Concatenate strings**

```r
cv2 <- paste(cv, 1:3, sep="")
```

Result:

```r
"a1" "b2" "c3"
```

### **Missing values**

```r
mvv <- c("a", "b", "c", NA)
is.na(mvv)
```

Result:

```r
FALSE FALSE FALSE TRUE
```

------

## **13. Type conversion**

```r
z <- 0:9
```

Creates:

```r
0 1 2 3 4 5 6 7 8 9
```

### **Check type**

```r
is.numeric(z)
```

Returns:

```r
TRUE
```

### **Convert to character**

```r
digits <- as.character(z)
```

Result:

```r
"0" "1" "2" ...
```

### **Convert back**

```r
d <- as.numeric(digits)
```

Result:

```r
0 1 2 3 ...
```

------

## **14. Subsetting vectors**

```r
x <- c(4, 5, 7, 3, 9)
```

### **Select positions 2 and 3**

```r
x[c(2,3)]
```

Result:

```r
5 7
```

### **Select values > 6**

```r
x[x > 6]
```

Result:

```r
7 9
```

### **Remove positions 1 and 5**

```r
x[-c(1,5)]
```

Result:

```r
5 7 3
```

------

## **15. Modify vectors**

```r
x[3] <- 23.1
```

Changes third value.

### **Append**

```r
c(x, 10.9)
```

### **Remove second element**

```r
x[-2]
```

### **Shorten vector**

```r
length(x) <- 2
```

Keeps only first two values.

------

## **16. Type coercion**

```r
c(1, "a")
```

Result:

```r
"1" "a"
```

Everything becomes character.

```r
c(1, TRUE)
```

Result:

```r
1 1
```

TRUE becomes 1.

```r
c(TRUE, "a")
```

Everything becomes character.

------

# **Data Frames**

## **17–19. Create and inspect data frames**

```r
df <- data.frame(
  name=c("Josh","rose"),
  age=c(23,35)
)
```

Result:

| **name** | **age** |
| -------- | ------- |
| Josh     | 23      |
| rose     | 35      |

### **First row**

```r
head(df,1)
```

### **Last row**

```r
tail(df,1)
```

### **Structure**

```r
str(df)
```

Output shows:

- name = character
- age = numeric

### **Access data**

```r
df[2,1]
```

Row 2, column 1.

```r
df[,2]
```

Entire second column.

------

# **Import and Export Data**

## **20. Read a tab-delimited file**

```r
cpm="https://.../cs.txt"
d <- read.delim(cpm)
head(d,3)
```

Downloads a text table from the web and shows first three rows.

Common RNA-seq workflow.

------

## **21. Write a file**

```r
x <- data.frame(a = "pi", b = pi)

write.table(
  x,
  file="foo.txt",
  sep="\t",
  row.names=FALSE
)
```

Creates:

```text
a     b
pi    3.141593
```

saved as `foo.txt`.

------

# **Graphics**

## **22. Scatter plot**

```r
area <- state.x77[, "Area"]
pop <- state.x77[, "Population"]
```

Built-in U.S. state dataset.

### **Plot**

```r
plot(area, pop)
```

Shows relationship between area and population.

### **Highlight Kansas**

```r
points(
 area["Kansas"],
 pop["Kansas"],
 col="purple",
 pch=19,
 cex=2
)
```

Adds a purple point.

------

## **23. Bar plot**

```r
barplot(
 pop/1000,
 las=2,
 cex.names=0.65
)
```

Arguments:

- `las=2` → vertical labels
- `cex.names` → label size

------

## **24. Histogram**

```r
hist(pop)
```

Shows distribution of state populations.

Useful for checking data distributions.

------

# **String Operations**

## **25. Count characters**

```r
cvec <- c("google", "hello", "the", "world")
nchar(cvec)
```

Result:

```r
6 5 3 5
```

------

## **26. Search strings**

```r
grep("o", cvec)
```

Finds strings containing `"o"`.



Returns positions:

```r
1 2 4
```

because

```text
google
hello
world
```

contain “o”.

------

## **27. Replace strings**

### **Replace first occurrence**

```r
sub("o", "O", cvec)
```

Result:

```r
"gOogle"
"hellO"
"the"
"wOrld"
```

Only first match is replaced.

### **Replace all occurrences**

```r
gsub("o", "O", cvec)
```

Result:

```r
"gOOgle"
"hellO"
"the"
"wOrld"
```

All matches are replaced.

------

## **28. Display variable**

```r
Sys.setenv("DISPLAY"=":0.0")
```

Sets the graphical display environment.

Mostly relevant when running R remotely on Linux systems with X11 graphics.

------

### **Key concepts covered in this tutorial**

1. Arithmetic operations
2. Variables and assignment
3. Vectors
4. Logical values
5. Character strings
6. Type conversion
7. Subsetting
8. Data frames
9. Reading and writing files
10. Basic plotting (scatter, barplot, histogram)
11. String searching and replacement

These are the foundational R skills needed before moving into RNA-seq analysis with packages such as **DESeq2** and **clusterProfiler**.
