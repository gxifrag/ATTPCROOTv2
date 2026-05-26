// Set a consistent, clean style
   { TStyle* mystyle = new TStyle("mystyle", "My Physics Style");

    // --- Canvas and Frame Settings ---
    mystyle->SetCanvasBorderMode(0);
    mystyle->SetPadBorderMode(0);
    mystyle->SetPadColor(0);
    mystyle->SetCanvasColor(0);
    mystyle->SetFrameFillColor(0);
    mystyle->SetFrameBorderMode(0);
    mystyle->SetPadTopMargin(0.12);    // Increase margin for title/stats text
    mystyle->SetPadBottomMargin(0.12); // Give space for X-axis title
    mystyle->SetPadLeftMargin(0.15);   // Give space for Y-axis title

    // --- Axis Settings ---
    mystyle->SetAxisColor(1, "XYZ");
    mystyle->SetLabelFont(42, "XYZ");
    mystyle->SetLabelSize(0.045, "XYZ");
    mystyle->SetTitleFont(42, "XYZ");
    mystyle->SetTitleFillColor(0);
    mystyle->SetTitleSize(0.04, "XYZ");
    mystyle->SetTitleOffset(1.3, "Y");   // Move Y-axis title away from label
    mystyle->SetTitleBorderSize(0);
    mystyle->SetTitleStyle(0);
    //mystyle->SetNdivisions(505, "XYZ");  // Remove ticks on top and right axes
    mystyle->SetTickLength(0.015, "XYZ");
    // --- Error Bar and Marker Settings ---
    mystyle->SetErrorX(0.0); // No horizontal caps on error bars
    mystyle->SetMarkerStyle(20); // Default marker style to filled circle (for TGraphErrors)

    mystyle->SetStatStyle(0); 
    // SetStatBorderSize(0) removes the box border.
    mystyle->SetStatBorderSize(0);
    // SetStatX/Y and W/H control position and size. Set them very small or outside the plot.
    mystyle->SetStatW(0.01); // Set width to minimum
    mystyle->SetStatH(0.01); // Set height to minimum

    // --- Legend and Stats Box Settings ---
    mystyle->SetLegendBorderSize(1);
    mystyle->SetStatBorderSize(0);
    mystyle->SetStatFont(42);
    mystyle->SetStatStyle(0); // Transparent background for stats box
    mystyle->SetOptStat(0);   // Globally disables drawing the stats box
    mystyle->SetStatFontSize(0.00); // Hides stats text if box somehow appears
    mystyle->SetLegendFont(42);
     mystyle->SetLegendFillColor(0); // Optionally set color to white/transparent
    mystyle->SetLegendTextSize(0.035); // NOTE: TStyle doesn't have SetLegendTextSize()
                                    // You must still set text size on the TLegend object itself.

    // --- Line and Fill Settings ---
    mystyle->SetLineWidth(2); // Thicker lines for curves

    gROOT->SetStyle("mystyle");
    gROOT->ForceStyle();
}
