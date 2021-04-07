dir = "/home/ula/catfish/test_scans/"

margin = 0;

scan_list = getFileList(dir);
for (n = 0; n < scan_list.length; n++){
	filename = scan_list[n];
	run("TIFF Virtual Stack...", "open=" + dir + filename);
	selectWindow(filename);

	noext = substring(filename, 0, filename.indexOf("."));

	splitDir = dir + noext +"/";
	File.makeDirectory(splitDir);

	slicenumber = nSlices/3;
	for (i = 0; i<3; i++) {
		start = i*slicenumber+1+margin;
		stop = (i+1)*slicenumber-margin;
		run("Make Substack...", "  slices="+start+"-"+stop);
		run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", splitDir + (i+1) + ".tif");
		close("Substack ("+start+"-"+stop + ")");
		selectWindow(filename);
	}
	close(filename);
	run("Images to Stack", "name=Stack title=[] use");
	run("Stack to RGB");
	saveAs("Tiff", splitDir + "orig.tif");
	run("Close All");
}



