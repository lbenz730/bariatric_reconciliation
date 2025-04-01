rsync -av -e ssh --exclude-from 'shell/exclude_list.txt' lbenz@fasselogin.rc.fas.harvard.edu:~/bariatric_tte/* .
Rscript scripts/copy_files.R
#rm Rplots.pdf