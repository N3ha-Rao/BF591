#BaseR

**Problem Statement**
This assignment will focus on basic functions of R with an emphasis on working with base R data structures and functional programming patterns.

**Learning Objectives and Skill List**
- Understand how vector, matrix, and data frame data structures operate in R
- Logical tests and operator functions applied to data structures
- Use functions and apply() to iterate over rows of a matrix
- Use functional programming concepts to generalize data summarizatoin operations
- Bonus: learn how to appropriately handle missing values

**Suggested workflow for developing and checking the code in this assignment:**
- main.R contains function definitions, including signature descriptions, for a number of functions, but the bodies of those functions are currently blank
- report.Rmd has code chunks that call functions defined in main.R - you do not need to write anything in the Rmd file (but you may)
- Your task is to read the function descriptions and the text in the Rmarkdown document and fill in the function bodies to produce the desired behavior in main.R
- You can test your work by executing individual code chunks in report.Rmd and comparing your output to the example compiled report in the repo
- In the workflow, you will go back and forth between developing code in main.R and running code chunks in report.Rmd
- In addition to inspecting your report results, also run testthat:test_file('test_main.R') to ensure they work correctly.
- When you have developed function bodies for all the functions and executed all the code chunks in the report successfully, you should be able to knit the entire report
