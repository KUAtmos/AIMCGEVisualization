#This program aims to make ppt file from selected indicators listed in ../data/pptlist.txt
# Packages
library(officer)
library(xml2)
# ----- Settings -----
txt_file <- "../data/pptlist.txt"      # file list with columns: name,left,top,width,height
pptx_out <- paste0(outdir,"IAMC.pptx")       # output pptx filename
layout_nm  <- "Title and Content" # slide layout
master_nm  <- "Office Theme"      # slide master
tpl_path   <- "../data/template_16x9.pptx"      # your 16:9 template saved from PowerPoint
# --------------------

# Helper: convert SVG length string to inches (supports in, cm, mm, px, pt, pc)
to_inches <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  x <- trimws(x)
  # Extract numeric value and unit
  m <- regexec("^([0-9.+-eE]+)\\s*([a-zA-Z%]*)$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) < 3) stop("Cannot parse SVG length: ", x)
  val <- as.numeric(r[2])
  unit <- tolower(r[3])
  
  # CSS/SVG defaults: px assumed when unit is missing
  if (unit %in% c("", "px"))      return(val / 96)        # 96 px per inch (CSS)
  if (unit == "in")               return(val)
  if (unit == "cm")               return(val / 2.54)
  if (unit == "mm")               return(val / 25.4)
  if (unit == "pt")               return(val / 72)        # 72 pt per inch
  if (unit == "pc")               return(val / 6)         # 1 pc = 12 pt
  stop("Unsupported SVG unit: ", unit)
}

# Helper: read nominal width/height (inches) from an SVG file
# - Prefer width/height attributes if present
# - Fallback to viewBox (pixels) assuming 96 px/inch
svg_nominal_size_in <- function(svg_path) {
  doc <- read_xml(svg_path)
  svg_node <- xml_find_first(doc, "/*[local-name()='svg']")
  if (is.na(svg_node)) stop("Not an SVG root: ", svg_path)
  
  w_attr <- xml_attr(svg_node, "width")
  h_attr <- xml_attr(svg_node, "height")
  
  if (!is.na(w_attr) && !is.na(h_attr) && nchar(w_attr) && nchar(h_attr)) {
    w_in <- to_inches(w_attr)
    h_in <- to_inches(h_attr)
    if (is.finite(w_in) && is.finite(h_in)) return(c(w_in, h_in))
  }
  
  # Fallback to viewBox: "minx miny width height"
  vb <- xml_attr(svg_node, "viewBox")
  if (!is.na(vb) && nchar(vb)) {
    parts <- as.numeric(strsplit(vb, "[ ,]+")[[1]])
    if (length(parts) == 4) {
      vw <- parts[3]; vh <- parts[4]
      return(c(vw, vh) / 96)  # assume 96 px/in
    }
  }
  
  stop("Cannot infer SVG size for: ", svg_path,
       " (missing width/height and viewBox)")
}

# Read file list (tab-style). If your file is whitespace-separated, use read.table(..., sep = "")
df <- read.table(txt_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Basic validation
required_cols <- c("name","title","regionflag","region", "left", "top", "scale","newslide")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) stop("../data/pptlist.txt: ", paste(missing_cols, collapse = ", "))

# Ensure numeric columns are numeric
num_cols <- c("left", "top", "scale","newslide")
df[num_cols] <- lapply(df[num_cols], as.numeric)

# Create document from 16:9 template
doc <- read_pptx(tpl_path)

# Loop over rows: each row -> one slide
for (i in seq_len(nrow(df))) {
  # Build SVG filepath from base name
  base   <- df$name[i]
  if(df$regionflag[i]=="Single"){
    svgf <- paste0(outdir,"byRegion/",df$region[i],"/svg/",df$name[i], "_",df$region[i],".svg")}
  if(df$regionflag[i]=="Multi"){
    svgf <- paste0(outdir,"multiReg",df$region[i],"/svg/",df$name[i], "_",df$region[i],".svg")}
  # Output pptx file
  pptx_out <- paste0(outdir,"iamc",Figureproj,".pptx")

  if (file.exists(svgf)) {
    # Extract geometry (inches)
    left   <- df$left[i]
    top    <- df$top[i]
    scale  <- df$scale[i]  
    # Keep aspect ratio by uniform scaling
    sz    <- svg_nominal_size_in(svgf)  # c(w,h) in inches
    w_out <- sz[1] * scale
    h_out <- sz[2] * scale
  
    if(df$newslide[i]==1){
      doc <- add_slide(doc, layout = layout_nm, master = master_nm)
      doc <- ph_with(doc, paste0(df$title[i]," ",df$region[i]), location = ph_location_type(type = "title"))
      #Add text box
      nbsp <- "\u00A0"  # non-breaking space
      doc <- ph_with(
        doc, value = nbsp,
        location = ph_location(left = 1, top = 5, width = 10, height = 2)
      )
    }
    doc <- ph_with(
      doc, external_img(svgf),
      location = ph_location(left = left, top = top, width = w_out, height = h_out)
    )
  }
}

# Save the pptx
print(doc, target = pptx_out)
