# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 19:02:37 2024

@author: andrea
"""

#Python script for running peak identification with matchms

from matchms.importing import load_from_msp
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import add_retention_index
from matchms.filtering import add_retention_time
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
import openpyxl

# Download msp file from MassBank of North America repository at https://mona.fiehnlab.ucdavis.edu/
file_msp_ref = "C:\\Users\\andre\\Desktop\\SUMMER_COURSE\\MONA\\MoNA-export-GC-MS_Spectra.msp"
references = list(load_from_msp(file_msp_ref))

file_msp = "C:\\Users\\andre\\Documents\\spectra\\msfinder_result.msp"
spectrums_raw = list(load_from_msp(file_msp))



spectrums = []
for spectrum in spectrums_raw:
    add_retention_time(spectrum)
    add_retention_index(spectrum)
    # Apply default filter to standardize ion mode, correct charge and more.
    # Default filter is fully explained at https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html .
    spectrum = default_filters(spectrum)
    # Scale peak intensities to maximum of 1
    spectrum = normalize_intensities(spectrum)
    spectrums.append(spectrum)
    

# Calculate Cosine similarity scores between all spectrums
# For other similarity score methods see https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html .
scores = calculate_scores(references=references,
                          queries=spectrums,
                          similarity_function=CosineGreedy(),
                          is_symmetric=False)

# This computed all-vs-all similarity scores, the array of which can be accessed as scores.scores
print(f"Size of matrix of computed similarities: {scores.scores.shape}")

# Matchms allows to get the best matches for any query using scores_by_query
query = spectrums[1]  # just an example
best_matches = scores.scores_by_query(query, 'CosineGreedy_score', sort=True)

# Print the calculated scores for each spectrum pair

for (reference, score) in best_matches[:1]:
    # Ignore scores between same spectrum
    if reference is not query:
        print(f"Reference compound name: {reference.metadata['compound_name']}")
        print(f"Score: {score[0]:.4f}")
        print(f"Number of matching peaks: {score[1]}")
        print("----------------------------")


#Get compound names
names = []
for i in range(len(spectrums)):
    query = spectrums[i]
    best_matches = scores.scores_by_query(query, 'CosineGreedy_score', sort=True)
    for (reference, score) in best_matches[:1]:
        if reference is not query:
            names.append(reference.metadata['compound_name'])








# to open the excel sheet and if it has macros
# srcfile = openpyxl.load_workbook('C:\\Users\\andre\\Documents\\datasets\\SpecAbund.xlsx', read_only=False, keep_vba=True)

# # get sheetname from the file
# sheetname = srcfile.get_sheet_by_name('SpecAbund')
# # write something in B2 cell of the supplied sheet
# sheetname['D1'] = names[0]
# # write to row 1,col 1 explicitly, this type of writing is useful to
# # write something in loops
# sheetname.cell(row=1, column=1).value = 'something'

# # save it as a new file, the original file is untouched and here I am saving
# # it as xlsm(m here denotes macros).
# srcfile.save('C:\\Users\\andre\\Documents\\datasets\\SpecAbund_matchms.xlsx')


    
