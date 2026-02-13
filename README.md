------------------------------------------------------------------------
FROM README Similarity/QTOF
------------------------------------------------------------------------

**The QTOF folder contains the following files :**

-   "sim_functions_qtofneg.R" : contains all the similarity functions.

-   "QTOF_neg_call.qmd" : in this file we call the similarity functions from the "sim_functions_qtofneg.R" to calculate the similarity between flax qtof negative and Dynlib qtof negative data.

-   "test_real_spectra.qmd" : here we test the similarity functions from "sim_functions_qtofneg.R" on some spectra examples and we compare the results with other similarity measures such as compareSpectra() from the spectra package, and NeutralLossesCosine() from the Python's matchms library.

-   "unit_test_qtofneg.R" : We used the "testthat" library to do a unit test for the functions in the "sim_functions_qtofneg.R" file.

-   "flax_qtof_to_Dynlib.qmd" : In this file we added the flax qtof negative data to the sql database of the Dynlib qtof negative.

**The execution order:**

1.  "sim_functions_qtofneg.R"

2.  "QTOF_neg_call.qmd"

3.  "unit_test_qtofneg.R"

4.  "test_real_spectra.qmd"

5.  "flax_qtof_to_Dynlib.qmd"

------------------------------------------------------------------------
FROM README Similarity/FTMS
------------------------------------------------------------------------
**The FTMS folder contains the following files :**

-   "sim_functions_ftmsneg.R" : contains all the similarity functions.

-   "FTMS_neg_call.qmd" : in this file we call the similarity functions from the "sim_functions_ftmsneg.R" to calculate the similarity between flax ftms negative and Dynlib ftms negative data.

-   "test_real_spectra.qmd" : here we test the similarity functions from "sim_functions_ftmsneg.R" on some spectra examples and we compare the results with other similarity measures such as compareSpectra() from the spectra package, and NeutralLossesCosine() from the Python's matchms library.

-   "unit_test_ftmsneg" : We used the "testthat" library to do a unit test for the functions in the "sim_functions_ftmsneg.R" file.

-   "flax_ftms_to_Dynlib.qmd" : In this file we added the flax ftms negative data to the sql database of the Dynlib ftms negative.

**The execution order:**

1.  "sim_functions_ftmsneg.R"

<!-- -->

1.  "FTMS_neg_call.qmd"

2.  "unit_test_ftmsneg"

3.  "test_real_spectra.qmd"

4.  "flax_ftms_to_Dynlib.qmd"

------------------------------------------------------------------------
# RDynLib

Creating and expanding an annotated spectral database generated from in house
LC-MS/MS and LC-MSn experiments.

In this repository we restructure the existing dynlib database that was on csv
format and we convert it to sql subdatabases in the following files:

- [dynlib-to-rforms-ftms-neg.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/dynlib-to-rforms-ftms-neg.qmd):
  in this file we convert the cvs format of the dynlib ftms negative database to
  a sql format

<!-- -->

- [dynlib-to-rforms-ftms-pos.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/dynlib-to-rforms-ftms-neg.qmd) :
  in this file we convert the cvs format of the dynlib ftms positive database to
  a sql format

- [dynlib-to-rforms-qtof-neg.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/dynlib-to-rforms-qtof-neg.qmd) :
  in this file we convert the cvs format of the dynlib qtof negative database to
  a sql format

- [dynlib-to-rforms-qtof-pos.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/dynlib-to-rforms-qtof-pos.qmd) :
  in this file we convert the cvs format of the dynlib qtof positive database to
  a sql format

- [DyL_SQL_to_MongoEmbadded2.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/DyL_SQL_to_MongoEmbadded2.qmd) :
  In this file, we created four MongoDB databases from the previous four SQL
  DynLib databases, along with a benchmark comparing the performance of both
  database types.

- [Similarity](https://github.com/rebdau/RDynLib/tree/dev_Ahlam/Similarity) :
  this folder contains two sub folders in which we calculate the similarity
  between flax data and dynlib database, in the first subfolder "FTMS" we
  calculate the similarity between flax ftms negative data and dynlib ftms
  negative data, in addition we created a script for adding ftmsneg flax to the
  dynlib ftmsneg sql database. The second folder "QTOF" we calculate the
  similarity between flax qtof negative data and dynlib qtof negative data, in
  addition we created a script for adding qtofneg flax to the dynlib qtofneg sql
  database.

 - [dyL_annot_FT_QT_pos.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/dyL_annot_FT_QT_pos.qmd) :
  This Quarto document represents the workflow of adding the annotation and grouping columns from
  the xcms2.txt files in the old dynlib structure to the new SQL version, specifically to the
  ms_compound table. In this file, we processed the dynlib FTMS positive data; the QTOF positive
  data does not contain any annotation or grouping file.

- [dyL_flax_annot_FTMS_neg.qmd](https://github.com/rebdau/RDynLib/blob/dev_Ahlam/dyL_flax_annot_FTMS_neg.qmd) :
  This Quarto document represents the workflow of adding the annotation and grouping columns from
  the xcms2.txt files in the old dynlib structure to the new SQL version, specifically to the
  ms_compound table. In this file, we processed the dynlib FTMS and QTOF negative data; Plus the annotation
  and grouping columns from the `XCMSExperiment` objects of both QTOF and FTMS negative flax data to the same
  sql database within Dynlib.
  

    
