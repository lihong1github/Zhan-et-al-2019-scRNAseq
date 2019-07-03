title = getTitle();
name1 = substring(title, 0, lastIndexOf(title, "."));

run("Input/Output...", "jpeg=100 gif=-1 file=.csv copy_column save_column");
//run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
run("Colors...", "foreground=white background=black selection=yellow");
m = getWidth;
n = getHeight;
dir = "/Users/lzhan/Desktop/Output/";
marker1 = "Iba1";
marker2 = "Mac2";
marker3 = "EdU";

selectWindow(title);
run("Crop");

run("Split Channels");

// identify Iba1 cells
selectWindow("C4-" + title);
run("Duplicate...", " ");
run("Smooth");
setMinAndMax(5, 150); //100
run("Apply LUT");
run("adaptiveThr ", "using=Mean from=100 then=-15");
//setAutoThreshold("MaxEntropy dark");
run("Convert to Mask");
run("Options...", "iterations=10 count=5 black do=Open");
run("Set Measurements...", "area redirect=None decimal=3");
run("Analyze Particles...", "size=50-500 show=Masks exclude add")
run("Convert to Mask");
run("Options...", "iterations=5 count=3 black do=Dilate");
rename("MARKER1_ROI");
run("Duplicate...", " ");
run("Outline");
rename("MARKER1_ROI_OUTLINE");
selectWindow("C4-" + title);
run("Duplicate...", " ");
run("Smooth");
setMinAndMax(5, 150); //100
run("Apply LUT");
rename("MARKER1");
run("Merge Channels...", "c1=MARKER1_ROI_OUTLINE c4=MARKER1 create keep ignore");
saveAs("Jpeg", dir + name1 + "_" + marker1);


// identify Mac2 cells
selectWindow("C2-" + title);
run("Duplicate...", " ");
run("Smooth");
setMinAndMax(5, 150); //100
run("Apply LUT");
run("adaptiveThr ", "using=Mean from=100 then=-20");
run("Convert to Mask");
run("Options...", "iterations=10 count=5 black do=Open");
run("Analyze Particles...", "size=50-800 show=Masks exclude");
run("Invert");
run("Convert to Mask");
rename("MARKER2_ROI");
run("Duplicate...", " ");
run("Options...", "iterations=5 count=3 black do=Dilate");
run("Outline");
rename("MARKER2_ROI_OUTLINE");
selectWindow("C2-" + title);
run("Duplicate...", " ");
run("Smooth");
setMinAndMax(5, 100); //100
run("Apply LUT");
rename("MARKER2");
run("Merge Channels...", "c2=MARKER2_ROI_OUTLINE c4=MARKER2 create keep ignore");
saveAs("Jpeg", dir + name1 + "_" + marker2);

// identify EdU cells
selectWindow("C3-" + title);
run("Duplicate...", " ");
run("Smooth");
setMinAndMax(0, 200); //100
run("Apply LUT");
run("adaptiveThr ", "using=Mean from=50 then=-20");
//setAutoThreshold("MaxEntropy dark");
run("Convert to Mask");
run("Set Measurements...", "area redirect=None decimal=3");
run("Analyze Particles...", "size=20-200 circularity=0.4-1 show=Masks exclude");
run("Invert");
run("Convert to Mask");
rename("MARKER3_ROI");
run("Duplicate...", " ");
run("Options...", "iterations=10 count=3 black do=Dilate");
run("Outline");
rename("MARKER3_ROI_OUTLINE");
selectWindow("C3-" + title);
run("Duplicate...", " ");
run("Smooth");
setMinAndMax(0, 200); 
run("Apply LUT");
rename("MARKER3");
run("Merge Channels...", "c1=MARKER3_ROI_OUTLINE c4=MARKER3 create keep ignore");
saveAs("Jpeg", dir + name1 + "_" + marker3);


run("Set Measurements...", "area_fraction redirect=None decimal=3");
run("Clear Results");

count=roiManager("count");

array_select = newArray();
array_all = newArray();
if(count==0) {
		run("Set Measurements...", "area area_fraction redirect=None decimal=3");
		run("Clear Results");
		saveAs("Results", dir + name1 + "_quantification.csv");
			list = getList("window.titles"); 
			for (i=0; i<list.length; i++){ 
			winame = list[i]; 
			selectWindow(winame); 
			run("Close"); } 
			run("Close All");
			}
else {

for (i=0;i<count;i++) {
	array_all = Array.concat(array_all,i);
	selectWindow("MARKER2_ROI"); 
	roiManager("Select", i);
	roiManager("Measure");
	thr=getResult("%Area", i);
	if(thr>10) {array_select = Array.concat(array_select,i);}
		}
	}

if(array_select.length==0) {
		run("Set Measurements...", "area area_fraction redirect=None decimal=3");
		run("Clear Results");
		saveAs("Results", dir + name1 + "_quantification.csv");
		
		run("Merge Channels...", " c1=MARKER2 c2=MARKER3 c4=MARKER1 create keep ignore");
 		run("RGB Color");
		saveAs("Jpeg", dir + name1 + "_" + marker1 + "_" + marker2 + "_" + marker3);
		
		list = getList("window.titles"); 
			for (i=0; i<list.length; i++){ 
			winame = list[i]; 
			selectWindow(winame); 
			run("Close"); } 
			run("Close All");
			}

	else {

// Save EdU+Iba1-ROI image for quality control

run("Set Measurements...", "area area_fraction redirect=None decimal=3");
run("Clear Results");
selectWindow("MARKER3_ROI");
roiManager("Select", array_select);
roiManager("Measure");
saveAs("Results", dir + name1 + "_quantification.csv");

run("Merge Channels...", " c1=MARKER2 c2=MARKER3 c4=MARKER1 create keep ignore");
run("RGB Color");
roiManager("Set Line Width", 1);
roiManager("select", array_select); 
roiManager("Draw");
saveAs("Jpeg", dir + name1 + "Iba1+Mac2+_select_EDU");

list = getList("window.titles"); 
			for (i=0; i<list.length; i++){ 
			winame = list[i]; 
			selectWindow(winame); 
			run("Close"); 
			} 
run("Close All");
	}




