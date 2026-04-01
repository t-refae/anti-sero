## this functionality is outside targets pipeline because 
## gt_save() cannot be used with gt_group objects for some reason

library(targets)
library(gt)
tar_source()

cva6_tab <- tar_read(CVA6_param_table)
ev71_tab <- tar_read(EV71_param_table)
ev68_tab <- tar_read(EV68_param_table)

combo_tab <- gt_group(cva6_tab, ev71_tab, ev68_tab)

combo_tab

# manually save pdf via html
