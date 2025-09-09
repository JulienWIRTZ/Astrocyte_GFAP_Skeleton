// ===== Macro for Analysis of astrocytes network [GFAP] (density, intensity, morphometry) =====
// Does not work for z-stack images, ONLY simple plan images. It is not necessary to define a ROI_background to substract background
// It is necessary to check the results of skeleton plan, a bad threshold can generate mistakes !!

requires("1.53");

// ------------------ Uilities functions ------------------

// Count 8-neighborhood degree on a binary image (0/255)
function _countNeighbors(title, x, y){
    selectWindow(title);
    n = 0;
    for (dy=-1; dy<=1; dy++){
        for (dx=-1; dx<=1; dx++){
            if (dx==0 && dy==0) continue;
            if (getPixel(x+dx, y+dy) > 0) n++;
        }
    }
    return n;
}

// Remove short spurs (branches that start at an endpoint and are <= Lmin_px long)
function pruneShortSpurs(srcTitle, Lmin_px, border_px){
    selectWindow(srcTitle); run("Duplicate...", "title=SKEL_work");
    selectWindow("SKEL_work");
    getDimensions(w,h,ch,z,t);

    for (iter=0; iter<20; iter++){
        changed = false;

        for (y=1; y<h-1; y++){
            for (x=1; x<w-1; x++){
                if (getPixel(x,y)==0) continue;

                // distance to edge
                dEdge = x;
                if (y < dEdge) dEdge = y;
                if ((w-1-x) < dEdge) dEdge = (w-1-x);
                if ((h-1-y) < dEdge) dEdge = (h-1-y);
                if (dEdge <= border_px) continue;

                // endpoint
                deg = 0;
                for (dy=-1; dy<=1; dy++) for (dx=-1; dx<=1; dx++){
                    if (dx==0 && dy==0) continue;
                    if (getPixel(x+dx, y+dy)>0) deg++;
                }
                if (deg != 1) continue;

                // ---- 1) measure length ----
                prevX = x; prevY = y; curX = x; curY = y; steps = 0;
                while (true){
                    // find unique neigh (≠ previous one)
                    neighCount = 0; nx = -1; ny = -1;
                    for (ddy=-1; ddy<=1; ddy++) for (ddx=-1; ddx<=1; ddx++){
                        if (ddx==0 && ddy==0) continue;
                        tx = curX+ddx; ty = curY+ddy;
                        if (getPixel(tx,ty)>0 && !(tx==prevX && ty==prevY)){
                            neighCount++; nx = tx; ny = ty;
                        }
                    }
                    if (neighCount==0) break;
                    if (neighCount>=2) break;
                    steps++;
                    if (steps > Lmin_px) break; // long enough -> keep it
                    prevX = curX; prevY = curY; curX = nx; curY = ny;
                }

                // ---- 2) if branche too short -> delete it ----
                if (steps <= Lmin_px){
                    prevX = x; prevY = y; curX = x; curY = y;
                    setPixel(curX, curY, 0);
                    while (true){
                        neighCount = 0; nx = -1; ny = -1;
                        for (ddy=-1; ddy<=1; ddy++) for (ddx=-1; ddx<=1; ddx++){
                            if (ddx==0 && ddy==0) continue;
                            tx = curX+ddx; ty = curY+ddy;
                            if (getPixel(tx,ty)>0 && !(tx==prevX && ty==prevY)){
                                neighCount++; nx = tx; ny = ty;
                            }
                        }
                        if (neighCount!=1) break;
                        setPixel(nx, ny, 0);
                        prevX = curX; prevY = curY; curX = nx; curY = ny;
                    }
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    return "SKEL_work";
}

// Count endpoints (excluding those within border_px from edges) and junctions (deg>=3)
function countEPandJN(title, border_px){
    selectWindow(title);
    getDimensions(w,h,ch,z,t);
    ep = 0; jn = 0;
    for (y=1; y<h-1; y++){
        for (x=1; x<w-1; x++){
            if (getPixel(x,y)==0) continue;

            // distance to border
            dEdge = x;
            if (y < dEdge) dEdge = y;
            if ((w - 1 - x) < dEdge) dEdge = (w - 1 - x);
            if ((h - 1 - y) < dEdge) dEdge = (h - 1 - y);

            deg = _countNeighbors(title, x, y);
            if (deg==1){
                if (dEdge > border_px) ep++;
            } else if (deg>=3){
                jn++;
            }
        }
    }
    return newArray(ep, jn);
}

// ------------------ Basic Macro ------------------

origTitle = getTitle();
selectWindow(origTitle);
rename("SRC");

// Calibration (µm/px) (check pixel calibration on the original picture .czi) (VERY IMPORTANT)
Dialog.create("Calibration (µm/px)");
Dialog.addNumber("Pixel width (µm/px)", 0.345);
Dialog.addNumber("Pixel height (µm/px)", 0.345);
Dialog.addString("Unit", "micron");
Dialog.show();
pw_in = Dialog.getNumber(); ph_in = Dialog.getNumber(); unit_in = Dialog.getString();
run("Properties...", "unit="+unit_in+" pixel_width="+pw_in+" pixel_height="+ph_in+" voxel_depth=1");
getPixelSize(unit, pw, ph, pd);
if (pw <= 0 || ph <= 0) exit("Invalid calibration.");

// Load ROIs
roiPath = File.openDialog("Choose ROIset (.zip)");
if (roiPath=="") exit("No ROIset.");
roiManager("Reset"); roiManager("Open", roiPath);
if (roiManager("count")==0) exit("ROIset is empty.");

// Select GFAP channel in the window
selectWindow("SRC");
titlesBefore = getList("image.titles");
run("Split Channels");
titlesAfter  = getList("image.titles");
// collection of new window titles
splits = newArray();
for (i=0; i<lengthOf(titlesAfter); i++){
    found = false;
    for (j=0; j<lengthOf(titlesBefore); j++) if (titlesAfter[i]==titlesBefore[j]) found = true;
    if (!found) splits = Array.concat(splits, titlesAfter[i]);
}
if (lengthOf(splits) > 0){
    def = splits[0];
    for (k=0; k<lengthOf(splits); k++) if (indexOf(splits[k], "(red)")>=0) def = splits[k];
    Dialog.create("Pick GFAP channel window");
    Dialog.addChoice("Window:", splits, def);
    Dialog.show();
    gfapWin = Dialog.getChoice();
    selectWindow(gfapWin); rename("GFAP_raw");
} else {
    rename("GFAP_raw");
}

// Preprocess (keep 16-bit for intensity and accuracy) (Adapt rolling and sigma for better results)
run("Duplicate...", "title=GFAP_pre16");
selectWindow("GFAP_pre16");
run("Subtract Background...", "rolling=50");
run("Gaussian Blur...", "sigma=0.8");

// 8-bit copy for segmentation
run("Duplicate...", "title=GFAP_pre8");
selectWindow("GFAP_pre8");
run("8-bit");

// Segmentation choice (other segmentations can be add in this list, Phansalkar can be adjust)
Dialog.create("Segmentation");
Dialog.addChoice("Method:", newArray("Otsu (global)", "Triangle (global)", "Phansalkar (local)"), "Otsu (global)");
Dialog.addNumber("Local radius (px) for Phansalkar", 15);
Dialog.show();
segChoice = Dialog.getChoice();
locRad    = Dialog.getNumber();

run("Duplicate...", "title=GFAP_mask");
selectWindow("GFAP_mask");
if (segChoice=="Phansalkar (local)"){
    run("Auto Local Threshold", "method=Phansalkar radius="+locRad+" parameter_1=0 parameter_2=0 white");
} else if (segChoice=="Triangle (global)"){
    setAutoThreshold("Triangle dark"); setOption("BlackBackground", true); run("Convert to Mask");
} else {
    setAutoThreshold("Otsu dark"); setOption("BlackBackground", true); run("Convert to Mask");
}

// Skeleton options (can be adjusted for better accuracy)
spurMin_um = 2.0;  // remove branches shorter than X µm
border_um  = 1.0;  // ignore endpoints close to ROI edge (distance is less than X µm)
pxSize = pw;       // assume square pixels
spurMin_px = round(spurMin_um / pxSize); if (spurMin_px < 1) spurMin_px = 1;
border_px  = round(border_um  / pxSize); if (border_px  < 0) border_px  = 0;

// Create output table
sumTab = "GFAP_DG_Summary_v15";
if (isOpen(sumTab)) Table.reset(sumTab); else Table.create(sumTab);

// Loop over ROIs
nR = roiManager("count");
for (ri=0; ri<nR; ri++){
    roiName = "ROI_" + ri;

    // %Area (mask)
    selectWindow("GFAP_mask");
    roiManager("Select", ri);
    getRawStatistics(nPixMask, meanMask, minMask, maxMask, stdMask);
    areaFrac = (meanMask/255.0) * 100.0;

    // Intensities (16-bit preprocessed)
    selectWindow("GFAP_pre16");
    roiManager("Select", ri);
    getRawStatistics(nPix, mean16, min16, max16, std16);
    intDen16 = mean16 * nPix;

    // ROI area
    area_um2 = nPix * pw * ph;

    // Skeletonize inside ROI (pruning & exclusion edge) ----
selectWindow("GFAP_mask"); run("Duplicate...", "title=SKEL_tmp");
selectWindow("SKEL_tmp");
roiManager("Select", ri); run("Clear Outside");
run("Skeletonize");

// 1) Pruning -> the function goes to "SKEL_work"
prunedTitle = pruneShortSpurs("SKEL_tmp", spurMin_px, border_px);

// 2) Normalize the name of Pruned window (depending on FIJI version)
if (!isOpen(prunedTitle)) prunedTitle = "SKEL_work";
if (!isOpen(prunedTitle)) prunedTitle = "SKEL_tmp";
if (!isOpen(prunedTitle)) exit("Aucune fenêtre de squelette trouvée après pruning.");
selectWindow(prunedTitle);
rename("SKEL_PRUNED");
print("Pruned skeleton window:", prunedTitle, "-> SKEL_PRUNED");

// 3) Length (µm)
getPixelSize(unit, pw, ph, pd); pxSize = pw;
selectWindow("SKEL_PRUNED");
getRawStatistics(nPixSkel, meanSkel, minSkel, maxSkel, stdSkel);
nonZero   = (meanSkel/255.0) * nPixSkel;
length_um = nonZero * pxSize;

// 4) Endpoints & Junctions (endpoints except on edge)
resEJ     = countEPandJN("SKEL_PRUNED", border_px);
endpoints = resEJ[0];
junctions = resEJ[1];

// Densities
area_mm2  = area_um2 / 1.0e6;
denom     = area_mm2; if (denom < 1.0e-12) denom = 1.0e-12;
lenDensity = (length_um / 1000.0) / denom;  // (mm/mm²) or mm-1
ep_per_mm2 = endpoints / denom;
jn_per_mm2 = junctions / denom;

// Mean thickness (µm) ≈ GFAP+ area (µm²) / length (µm)
area_GFAP_um2 = (areaFrac/100.0) * area_um2;
if (length_um > 0)
    mean_thick_um = area_GFAP_um2 / length_um;
else
    mean_thick_um = 0;

// Normalize on length (of mm of GFAP skeleton)
len_mm = length_um / 1000.0;
if (len_mm > 0) {
    tips_per_mm     = endpoints / len_mm;
    branches_per_mm = junctions / len_mm;
} else {
    tips_per_mm     = 0;
    branches_per_mm = 0;
}

// 6) Cleaning ?
close("SKEL_tmp");
close("SKEL_PRUNED");

    // QC overlay (PNG file)
    selectWindow("GFAP_pre16"); run("Duplicate...", "title=QC_base");
    selectWindow("QC_base"); run("Enhance Contrast", "saturated=0.35 normalize"); run("8-bit");
    selectWindow("GFAP_mask"); run("Duplicate...", "title=QC_mask");
    selectWindow("QC_mask"); roiManager("Select", ri); run("Clear Outside");
    run("RGB Merge...", "red=QC_mask green=QC_base");
    rename("QC_merged");
    setForegroundColor(255,255,0);
    roiManager("Select", ri); run("Line Width...", "line=2"); run("Draw");
    close("QC_merged"); close("QC_base"); close("QC_mask");
    close("SKEL_tmp"); close("SKEL_PRUNED");

    // DATA SET ! The best :)
    row = Table.size(sumTab);
    Table.set("Image", row, origTitle);
    Table.set("ROI", row, roiName);
    Table.set("PixelWidth_um",  row, pw);
    Table.set("PixelHeight_um", row, ph);
    Table.set("Unit",           row, unit);
    Table.set("Area_um2",       row, area_um2);
    Table.set("GFAP_%Area",     row, areaFrac);
    Table.set("GFAP_Mean16",    row, mean16);
    Table.set("GFAP_IntDen16",  row, intDen16);
    Table.set("Skel_Length_um", row, length_um);
    Table.set("Skel_LengthDensity_mm_per_mm2", row, lenDensity);
    Table.set("Skel_Endpoints", row, endpoints);
    Table.set("Skel_Endpoints_per_mm2", row, ep_per_mm2);
    Table.set("Skel_Junctions", row, junctions);
    Table.set("Skel_Junctions_per_mm2", row, jn_per_mm2);
    Table.set("Skel_Tips_per_mm", row, tips_per_mm);
    Table.set("Skel_Branches_per_mm", row, branches_per_mm);
    Table.set("GFAP_MeanThickness_um", row, mean_thick_um);
}
Table.update(sumTab);
