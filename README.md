# epistasis

pipeline.py:
1) Open table with mutations and read mutations one by one.
2) Make directory to save results. Everything is saved in "result.tsv". For each mutation save pdb and file with energies from FoldX.
3) Move to /tmp not to care about unnecessary files created by FoldX and Eris.
4) Launch FoldX, which also checks if given AA really is in the given position, save dG and ddG in result.tsv
5) Launch Eris with timeout of 5 minutes.
