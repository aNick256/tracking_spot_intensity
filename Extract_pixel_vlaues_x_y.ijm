// Set the inner and outer radii
r = 4 ;
inner_r = 8;
outer_r = 9;
run("Select None");
getPixelSize(unit, pixelWidth, pixelHeight);
pixelWidth = pixelWidth * 1000;

signal_area = PI * r * r ;
outer_area = PI * outer_r * outer_r ;
inner_area = PI * inner_r * inner_r ;
background_area = outer_area - inner_area ;

getDimensions(width, height, channels, slices, frames);

n_slice = slices / channels ;
directory = getDirectory("image");
// Get the title of the image
imgID = getImageID();
img_t = getTitle();
img_tile_array = split(img_t , ".");
img_title = img_tile_array[0];

// Define the subfolder name
subfolder_name = img_title + "_trace";
if(channels>1){
run("Split Channels");
selectImage("C1-" + img_t);
close();}

//run("Enhance Contrast", "saturated=0.35");
img_t = getTitle();
img_tile_array = split(img_t , ".");
img_title = img_tile_array[0];
// Get the directory of the open image
//

new_Directory = directory + File.separator + subfolder_name + File.separator ;
// Create the subfolder if it doesn't exist
File.makeDirectory(new_Directory);
Coordinates_file = new_Directory + img_title + "_Coordinates.csv" ;
result_file_dir = new_Directory + img_title + "_Fiji_Resuls_Table.csv";
f = File.open(Coordinates_file);


// Set the directory to the subfolder
// setDirectory(directory + subfolder_name);
//run("Subtract Background...", "rolling=50 stack");
//run("Enhance Contrast", "saturated=0.35");
run("Gaussian Fit");
//run("Simple Fit", "use_current_calibration");
//run("Gaussian Fit", "smoothing=0.32 box_size=1 background=1600 min_height=300 fraction_above_background=0.02 min_width=4 top_n=0 block_find_algorithm border=0 psf=[Circular Gaussian 2D] fit_background max_iterations=20 relative_threshold=1.000E-5 absolute_threshold=1.000E-10 single_fit single_region_size=10 initial_stddev=0.000 log_progress filter_results show_fit");
xvalues = newArray();
yvalues = newArray();

selectWindow("Fit Results");
saveAs("Text", result_file_dir);
close("Fit Results");

     lineseparator = "\n";
     cellseparator = ",";

     lines=split(File.openAsString(result_file_dir), lineseparator);

     // recreates the columns headers
     labels=split(lines[0], cellseparator);
     if (labels[0]==" ")
        k=1; // it is an ImageJ Results table, skip first column
     else
        k=0; // it is not a Results table, load all columns
     for (j=k; j<labels.length; j++)
        setResult(labels[j],0,0);

     // dispatches the data into the new RT
     run("Clear Results");
     for (i=1; i<lines.length; i++) {
        items=split(lines[i], cellseparator);
        for (j=k; j<items.length; j++)
           setResult(labels[j],i-1,items[j]);
     }
     updateResults();

N_results = nResults();


// Loop through each row in the ResultsTable
for(i = 0; i < N_results; i++) {
    // Get X and Y values from columns
    xValue = getResult("X (px)", i);
    yValue = getResult("Y (px)", i);

    // Convert values to double and nm
    xValue_nm = parseFloat(xValue) ;
    yValue_nm = parseFloat(yValue)  ;
	
    // Store values in arrays
    xvalues[i] = xValue_nm * pixelWidth ;
    yvalues[i] = yValue_nm * pixelWidth ;
}
header = "# {'pix_x': " + width + ", 'pix_y': " + height + ", 'pix_size': " + pixelWidth + ", 'r_peak': " + r + ", 'r_bg1': " + inner_r + ", 'r_bg2': " + outer_r + ", 'min_dist': " + outer_r + ", 'binning': 1}" ;
            // Create the CSV content
            print(f, header);
            Content = " id,x [nm],y [nm], center_pix";
            for (k = 0; k < frames; k++) {
                Content = Content + "," + k;
            }
            print(f, Content);



// Loop through each row in the ResultsTable
for (i = 0; i < N_results; i++) {
    ttl = img_t;

    // Get X and Y values from arrays
    xValue_nm = xvalues[i] / pixelWidth;
    yValue_nm = yvalues[i] / pixelWidth;

    // Check if the circular ROI for signal falls within the image boundaries
    if (xValue_nm - outer_r >= 0 && yValue_nm - outer_r >= 0 && xValue_nm + outer_r < getWidth() && yValue_nm + outer_r < getHeight()) {

        // Check if the distance from other dots is greater than or equal to 2 * outer_r
        excludeSpot = false;
        for (j = 0; j < N_results; j++) {
            if (i != j) {
                // Calculate the distance between the current spot and other dots
                distance = Math.sqrt(Math.pow(xValue_nm - xvalues[j] / pixelWidth, 2) + Math.pow(yValue_nm - yvalues[j] / pixelWidth, 2));
                if (distance < 2 * outer_r) {
                    excludeSpot = true;
                    break;  // No need to check further, exclude the spot
                }
            }
        }

        if (!excludeSpot) {
            // Define the circular ROI for signal (radius 2)
            makeOval(xValue_nm - r, yValue_nm - r, 2 * r, 2 * r);
            roiManager("add");
            selectWindow(ttl);

            // Run "Measure" to obtain measurements within the circular ROI
            roiManager("Select", 0);
            run("Measure Stack...");

            n_of_frames = nResults();
            // Get the mean intensity value from the "Results" window
            signalIntensity = newArray();
            for (k = 0; k < n_of_frames; k++) {
                signalIntensity[k] = getResult("Max", k);
            }
            run("Clear Results");

            // Define the annular ROI for background (between inner_r and outer_r)
            makeOval(xValue_nm - outer_r, yValue_nm - outer_r, 2 * outer_r, 2 * outer_r);
            roiManager("add");
            selectWindow(ttl);
            roiManager("Select", 1);
            run("Measure Stack...");
            BackgroundIntensity2 = newArray();
            for (k = 0; k < n_of_frames; k++) {
                BackgroundIntensity2[k] = getResult("Mean", k);
            }
            run("Clear Results");

            selectWindow(ttl);
            makeOval(xValue_nm - inner_r, yValue_nm - inner_r, 2 * inner_r, 2 * inner_r);
            roiManager("add");
            // Run "Measure" again to obtain measurements within the annular ROI
            selectWindow(ttl);
            roiManager("Select", 2);
            run("Measure Stack...");
            BackgroundIntensity1 = newArray();
            for (k = 0; k < n_of_frames; k++) {
                BackgroundIntensity1[k] = getResult("Mean", k);
            }
            //close("Results");
            background = newArray(n_of_frames);
            correctedIntensity = newArray();
            min_intensity = 0;
            for (k = 0; k < n_of_frames; k++) {
                background[k] = (BackgroundIntensity2[k]*outer_area - BackgroundIntensity1[k]*inner_area)/background_area;

                // Subtract background from the signal intensity
                correctedIntensity[k] = signalIntensity[k] - background[k];
                
                if(correctedIntensity[k]<min_intensity){
                	min_intensity = correctedIntensity[k] ;
                }

                // Print or store the corrected intensity
            }
            
            for (k = 0; k < n_of_frames; k++) {
            	correctedIntensity[k] = correctedIntensity[k] - (min_intensity) ;
            	
            }
            print("signal int:" + signalIntensity[0] + "  corrected_signal:" + correctedIntensity[0] + "  backg1:" + BackgroundIntensity1[0] + "  backg2:" + BackgroundIntensity2[0] + "  background:" + background[0]) ; 
            close("Results");
            roiManager("reset");
            

			center = width * floor(yvalues[i] / pixelWidth) + floor(xvalues[i] / pixelWidth) ;

            Content = " " + i + "," + xvalues[i] + "," + yvalues[i] + "," + center;
            for (j = 0; j < n_of_frames; j++) {
                Content = Content + "," + signalIntensity[j];
            }
            print(f, Content);
        }
    }
}



print("Corrected Intensity at");


File.close(f);
close("Gaussian Fit");
close("Results");

close();