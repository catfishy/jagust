rm(list= ls())
library(oro.nifti)
library(tractor.base)
library(XLConnect)
library(gdata)
library(fields)

sids = c('B06-203','B06-235','B07-227','B07-238','B07-243','B07-276','B07-277','B08-221','B08-223','B09-211','B09-216','B09-277','B10-282','B10-298','B10-306','B10-311','B10-313','B11-211','B11-232','B11-243','B11-253','B11-255','B11-256','B11-274','B12-202','B12-203','B12-209','B12-219','B12-226','B12-276','B12-277','B12-305','B12-316','B13-204','B13-205','B13-207','B13-225','B13-227','B13-238','B13-267','B13-280','B13-288','B13-289','B13-307','B13-313','B13-314','B13-317','B13-318','B13-319','B13-321','B13-341','B13-349','B13-352','B13-357','B13-362','B14-224','B14-237','B14-243','B14-257','B14-265','B14-274','B14-284','B14-285','B14-289','B14-292','B14-301','B14-303','B14-304','B15-214','B15-218','B15-228','B15-230','B15-234','B15-235','B15-249','B15-251','B15-254','B15-257','B15-263','B15-264','B15-278','B15-281','B15-283','B15-284','B15-288','B15-291','B15-293','B15-295','B15-299','B15-303','B15-309','B15-315','B15-316','B15-319','B15-323','B15-324','B15-325','B15-327','B15-331','B16-205','B16-212','B16-213','B16-220','B16-226','B16-227','B16-229','B16-231','B16-232','B16-236','B16-240','B16-241','B16-242','B16-246','B16-247','B16-254','B16-262','B16-264','B16-265','B16-269','B16-270','B16-271','B16-272','B16-273','B16-274','B16-275','B16-278','B16-282','B16-283','B16-298')

# wd = "/home/jagust/adni/adni_av45/qc_output/"
# wd = "/home/jagust/adni/adni_av1451/qc_output/"
wd = "/home/jagust/ahorng/bacs_qc_output/"

setwd(wd)
getwd()

# INPUT PARAMS
aparc_roi = c('8','47'); # cereb gray
# aparc_roi = c('8','47','7','46'); # whole cerebellum
number_of_slice_display = 24; #prefer square values -1 such as 15 24 35
# nu_loc = '/mri/nu.nii'
# aparc_loc = '/mri/aparc+aseg.nii'
nu_loc = '/mri/rnu.nii'
aparc_loc = '/mri/raparc+aseg.nii'

# DATA INPUT - EXCEL
data_table = read.csv(paste0(wd,'qc_input.csv'))

# Filter by SID
data_table = data_table[data_table$RID %in% sids,]

nb_subjects = length(data_table$RID)
data_table$pet_path = paste0(data_table$pet_path)
data_table$fs_path = paste0(data_table$fs_path)
data_table$nu_path = paste0(data_table$fs_path,nu_loc)
data_table$aparc_path = paste0(data_table$fs_path,aparc_loc)

pdf_allsubj_fs = paste0("pdfjam -q --papersize '{5in,5in}' ")
pdf_allsubj_aparc_sagittal = paste0("pdfjam -q --papersize '{5in,5in}' ")
pdf_allsubj_coreg_axial = paste0("pdfjam -q --papersize '{5in,5in}' ")
pdf_allsubj_coreg_sagittal = paste0("pdfjam -q --papersize '{5in,5in}' ")

display_img = c(0.2,1)
max_display = max(display_img)
palette_display = designer.colors(n=256, col= c("red", "purple4", "#05004B", "blue3", "blue", "#105EEE", "#34B8AD", "#00FFA6", "#00FF37", "#00FF51", "#FFFF00", "#FFB300", "#FF5100", "#FF2200", "#EE1B1B", "#C91F1F", "#AE1717"),
                                          x=seq(0,1,1/16),alpha=0.5)

for (i in 1:nb_subjects) {   # loop all subjects
  name_pdf_subj = paste0(data_table$RID[i],'_','qc', ".pdf")
  
  # test the existence of the pdf
  if (!file.exists(name_pdf_subj)){
    
    print(paste0('Create PDF for subject ',i,'/',nb_subjects,' - ',data_table$lbl_id[i],' - ', Sys.time()))
    print(data_table$aparc_path[i])
    print(data_table$nu_path[i])
    print(data_table$pet_path[i])
    
	  ################## CEREB GRAY ON MRI / SAGITTAL
    if (file.exists(data_table$aparc_path[i]) && file.exists(data_table$nu_path[i])){
  		# GREP COORDINATES of the ROI in native space FOR THE DISPLAY
  		ROI_load_int = readNIfTI(data_table$aparc_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
  		ROI_load_int[is.na(ROI_load_int)] = 0
  		ROI_load = ROI_load_int
  		ROI_load[ROI_load_int %in% aparc_roi] = 1
  		ROI_load[!(ROI_load_int %in% aparc_roi)] = 0

  		Xleft = 1
  		while(sum(ROI_load[Xleft,,],na.rm = T) <= 0) {Xleft=Xleft+1}
  		Xright = ROI_load@dim_[2]
  		while(sum(ROI_load[Xright,,],na.rm = T) <= 0) {Xright=Xright-1}
  		Xint = round((Xright - Xleft)/number_of_slice_display,1)
  		XseqROI = seq(Xleft,Xright,Xint)
  		if(length(XseqROI)>(number_of_slice_display+1)){
  		  while(length(XseqROI) > (number_of_slice_display+1)){; Xint = Xint + 0.1;  XseqROI = seq(Xleft,Xright,Xint);}
  		}
  		if(length(XseqROI)==number_of_slice_display){XseqROI=c(XseqROI,Xright)}
  		
  		# Overlay of the ROI on native space MRI
  		ROI_load_display = ROI_load
  		ROI_load_display[ROI_load_display>max_display] = max_display  #step necessary only when values in the image are above the one displayed.
  		tiff("01.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
  		T1_load = readNIfTI(data_table$nu_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
  		overlay.nifti(T1_load,ROI_load_display, z = XseqROI, col.x = gray(0:64/64),col.y= palette_display, zlim.x = NULL, 
  			  zlim.y  = display_img, plane = "sagittal", plot.type = "single")
  		title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2,   font.main= 2, col.main= "White", line = 3)
  		title(sub = paste0(data_table$aparc_path[i],' - size ROI: ',length(ROI_load>0)), cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
  		dev.off()
    } else {
  		tiff("01.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
  		plot.new()
  		title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
  		title(sub = 'NO FS output for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
  		dev.off()
    }
	
	  ################## APARC REGIONS ON PET / SAGITTAL
    if (file.exists(data_table$aparc_path[i]) && file.exists(data_table$pet_path[i])){
      # GREP COORDINATES of the ROI in native space FOR THE DISPLAY
      ROI_load_int = readNIfTI(data_table$aparc_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
      ROI_load_int[is.na(ROI_load_int)] = 0
      ROI_load = ROI_load_int
      
      Xleft = 1
      while(sum(ROI_load[Xleft,,],na.rm = T) <= 0) {Xleft=Xleft+1}
      Xright = ROI_load@dim_[2]
      while(sum(ROI_load[Xright,,],na.rm = T) <= 0) {Xright=Xright-1}
      Xint = round((Xright - Xleft)/number_of_slice_display,1)
      XseqROI = seq(Xleft,Xright,Xint)
      if(length(XseqROI)>(number_of_slice_display+1)){
        while(length(XseqROI) > (number_of_slice_display+1)){; Xint = Xint + 0.1;  XseqROI = seq(Xleft,Xright,Xint);}
      }
      if(length(XseqROI)==number_of_slice_display){XseqROI=c(XseqROI,Xright)}
      
      # Overlay of the ROI on native space MRI
      ROI_load_display = ROI_load
      tiff("02.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      T1_load = readNIfTI(data_table$pet_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
      overlay.nifti(T1_load,ROI_load_display, z = XseqROI, col.x = gray(0:64/64),col.y= palette_display, zlim.x = NULL, 
                    zlim.y  = c(0.2,2040), plane = "sagittal", plot.type = "single")
      title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2,   font.main= 2, col.main= "White", line = 3)
      title(sub = paste0(data_table$aparc_path[i]), cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    } else if (!file.exists(data_table$aparc_path[i])) {
      tiff("02.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      plot.new()
      title(paste0(data_table$RID[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
      title(sub = 'NO APARC for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    } else {
      tiff("02.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      plot.new()
      title(paste0(data_table$RID[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
      title(sub = 'NO PET for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    }
    
    ################## PET ON MRI / AXIAL
    if (file.exists(data_table$nu_path[i]) && file.exists(data_table$pet_path[i])){
      # GREP COORDINATES of the ROI in native space FOR THE DISPLAY
      ROI_load_int = readNIfTI(data_table$pet_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
      ROI_load_int[is.na(ROI_load_int)] = 0
      ROI_load = ROI_load_int
      
      Zleft = 1
      while(sum(ROI_load[,,Zleft],na.rm = T) <= 0) {Zleft=Zleft+1}
      Zright = ROI_load@dim_[4]
      while(sum(ROI_load[,,Zright],na.rm = T) <= 0) {Zright=Zright-1}
      Zint = round((Zright - Zleft)/number_of_slice_display,1)
      ZseqROI = seq(Zleft,Zright,Zint)
      if(length(ZseqROI)>(number_of_slice_display+1)){
        while(length(ZseqROI) > (number_of_slice_display+1)){; Zint = Zint + 0.1;  ZseqROI = seq(Zleft,Zright,Zint);}
      }
      if(length(ZseqROI)==number_of_slice_display){ZseqROI=c(ZseqROI,Zright)}

      # Only the ROI by itself
      ROI_load_display = ROI_load
      tiff("03.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      image(ROI_load_display, z = ZseqROI, col = palette_display, zlim  = c(0.5,1.8), plane = "axial", plot.type = "single")
      title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2, font.main= 2, col.main= "White", line = 3)
      title(sub = paste0(data_table$pet_path[i]), cex.sub = 0.5, font.sub= 2, col.sub= "RED")
      dev.off()
      
      # Overlay of the ROI on native space MRI
      ROI_load_display = ROI_load
      tiff("04.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      T1_load = readNIfTI(data_table$nu_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
      overlay.nifti(T1_load,ROI_load_display, z = ZseqROI, col.x = gray(0:64/64),col.y= palette_display, zlim.x = NULL, 
                    zlim.y  = c(0.5,1.8), plane = "axial", plot.type = "single")
      title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2,   font.main= 2, col.main= "White", line = 3)
      title(sub = paste0(data_table$pet_path[i]), cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    } else {
      tiff("03.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      plot.new()
      title(paste0(data_table$RID[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
      title(sub = 'NO NU/PET for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
      tiff("04.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      plot.new()
      title(paste0(data_table$RID[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
      title(sub = 'NO NU/PET for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    }
    
    ################## PET ON MRI / SAGITTAL
    if (file.exists(data_table$nu_path[i]) && file.exists(data_table$pet_path[i])){
      # GREP COORDINATES of the ROI in native space FOR THE DISPLAY
      ROI_load_int = readNIfTI(data_table$pet_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
      ROI_load_int[is.na(ROI_load_int)] = 0
      ROI_load = ROI_load_int
      
      Zleft = 1
      while(sum(ROI_load[Zleft,,],na.rm = T) <= 0) {Zleft=Zleft+1}
      Zright = ROI_load@dim_[2]
      while(sum(ROI_load[Zright,,],na.rm = T) <= 0) {Zright=Zright-1}
      Zint = round((Zright - Zleft)/number_of_slice_display,1)
      ZseqROI = seq(Zleft,Zright,Zint)
      if(length(ZseqROI)>(number_of_slice_display+1)){
        while(length(ZseqROI) > (number_of_slice_display+1)){; Zint = Zint + 0.1;  ZseqROI = seq(Zleft,Zright,Zint);}
      }
      if(length(ZseqROI)==number_of_slice_display){ZseqROI=c(ZseqROI,Zright)}
      
      # Only the ROI by itself
      ROI_load_display = ROI_load
      tiff("05.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      image(ROI_load_display, z = ZseqROI, col = palette_display, zlim  = c(0.5,1.8), 
            plane = "sagittal", plot.type = "single")
      title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2, font.main= 2, col.main= "White", line = 3)
      title(sub = paste0(data_table$pet_path[i]), cex.sub = 0.5, font.sub= 2, col.sub= "RED")
      dev.off()
      
      # Overlay of the ROI on native space MRI
      ROI_load_display = ROI_load
      tiff("06.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      T1_load = readNIfTI(data_table$nu_path[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL)
      overlay.nifti(T1_load,ROI_load_display, z = ZseqROI, col.x = gray(0:64/64),col.y= palette_display, zlim.x = NULL, 
                    zlim.y  = c(0.5,1.8), plane = "sagittal", plot.type = "single")
      title(paste0(data_table$RID[i],'-',data_table$name[i]), cex.main = 1.2,   font.main= 2, col.main= "White", line = 3)
      title(sub = paste0(data_table$pet_path[i]), cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    } else {
      tiff("05.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      plot.new()
      title(paste0(data_table$RID[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
      title(sub = 'NO NU/PET for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
      tiff("06.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'none')
      plot.new()
      title(paste0(data_table$RID[i]), cex.main = 1.2,   font.main= 2, col.main= "Black", line = 3)
      title(sub = 'NO NU/PET for the subj', cex.sub = 0.5,   font.sub= 2, col.sub= "RED")
      dev.off()
    }
    
    
    ############################################
    # Create PDF subject #######################
    ############################################
	  #intermediary step using tiff to prevent very heavy pdf (scrolling works better that way)
    #create one tiff with all
    one_tif = paste0('tiffcp 01.tiff 02.tiff 03.tiff 04.tiff 05.tiff 06.tiff temp_all.tiff')
    system(one_tif)
    #create one pdf with all
    one_pdf_subj = paste0('tiff2pdf -z -o ',name_pdf_subj," temp_all.tiff")
    system(one_pdf_subj)
    #remove temp
    del_tem = paste0('rm -rf 01.tiff 02.tiff 03.tiff 04.tiff 05.tiff 06.tiff temp_all.tiff')
    system(del_tem)
    
  }  else {print(paste0('PDF already exist for subject ',i,'/',nb_subjects,' - ',data_table$lbl_id[i],' - ',Sys.time()))}
  
  pdf_allsubj_fs = paste0(pdf_allsubj_fs,name_pdf_subj, " 1 ")  
  pdf_allsubj_aparc_sagittal = paste0(pdf_allsubj_aparc_sagittal,name_pdf_subj, " 2 ")
  pdf_allsubj_coreg_axial = paste0(pdf_allsubj_coreg_axial,name_pdf_subj, " 4 ")
  pdf_allsubj_coreg_axial = paste0(pdf_allsubj_coreg_axial,name_pdf_subj, " 3 ")
  pdf_allsubj_coreg_sagittal = paste0(pdf_allsubj_coreg_sagittal,name_pdf_subj, " 6 ")
  pdf_allsubj_coreg_sagittal = paste0(pdf_allsubj_coreg_sagittal,name_pdf_subj, " 5 ")
}

pdf_allsubj_fs = paste0(pdf_allsubj_fs,"-o allsubj_fs.pdf")
pdf_allsubj_aparc_sagittal = paste0(pdf_allsubj_aparc_sagittal,"-o allsubj_aparc_sagittal.pdf")
pdf_allsubj_coreg_axial = paste0(pdf_allsubj_coreg_axial,"-o allsubj_coreg_axial.pdf")
pdf_allsubj_coreg_sagittal = paste0(pdf_allsubj_coreg_sagittal,"-o allsubj_coreg_sagittal.pdf")

system(pdf_allsubj_fs)
system(pdf_allsubj_aparc_sagittal)
system(pdf_allsubj_coreg_axial)
system(pdf_allsubj_coreg_sagittal)
