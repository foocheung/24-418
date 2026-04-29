## Notes and Corrections

### Figure labeling issue (Figure 5D–G)
The labels are currently incorrect—**β-glucan (BG)** and **SYKi** appear to be switched and should be corrected.

### Missing clarification in text
The manuscript does not state that only the **BG condition** yielded a sufficient number of significant peaks for downstream analysis (see below I spent some time making sure this is clear in the methods)

### Code availability
Add in the link to the available analysis code that was used in the MS for the scRNA and scATAC:  
https://github.com/foocheung/24-418  

### Methods documentation
The methods section has been rewritten and expanded, including the addition of **20 references** covering the code, tools, and analytical frameworks used in this study.  

The **BG condition is now explicitly described**, including how it was used in the peak-level **FDR correction and downstream filtering strategy**.  

Additionally, for **peak-to-gene (peak2gene) analyses**, the methods now clarify that:
- Significant peak-level signals were primarily observed in the **BG condition**,  
- Other training conditions (SYKi, MDP) had limited peaks passing stringent FDR thresholds,  
- Therefore, a **peak-to-gene linkage–based approach** was used to enable cross-condition comparisons by leveraging correlated accessibility rather than relying strictly on peak-level significance. 

Full methods write-up:  
https://github.com/foocheung/24-418/blob/main/Methods_v01.docx
