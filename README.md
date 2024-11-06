# **BF591: R for Biological Sciences**

Welcome to the **BF591: R for Biological Sciences** course repository! This repository contains all the assignments for the course, and each assignment follows a similar format and workflow. Below is a description of the assignment structure and how the files interact.

[Course overview](https://bu-bioinfo.github.io/r-for-biological-sciences/index.html)

---

### **File Structure Overview**

![image](https://github.com/user-attachments/assets/65b6f4b9-c8c8-48af-9af0-634214fc05ba)

The repository for each assignment will have the following files:

```
├── reference_report.html
├── main.R
├── README.md
├── report.Rmd
└── test_main.R
```

---

### **Description of Files**

- **`reference_report.html`**  
  This file contains the completed **report** for the assignment. The student's task will be to replicate this report, with some flexibility in elements like plots, captions, and other details.

- **`main.R`**  
  This is the primary R script for the assignment. It contains **function definitions and descriptions** that are initially empty. The student will need to complete these functions as described in the script and ensure that they pass all the provided tests.

- **`README.md`**  
  Important information about the repository, installation instructions, and helpful links back to the assignment details will be found here.
  
- **`report.Rmd`**  
  This is an **R Markdown** file where the student will source their `main.R` script. It is set up with empty code chunks, and the goal is to use the functions from `main.R` to generate the necessary plots, figures, and analysis to match the report in `reference_report.html`.

- **`test_main.R`**  
  This file contains a **testing script** that ensures the code is running as expected. It sources the `main.R` script, executing it and checking the functions and outputs. The tests will help identify any issues with your code and verify the results.

