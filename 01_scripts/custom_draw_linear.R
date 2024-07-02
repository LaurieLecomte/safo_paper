# A custom implementation of the draw_linear() function originally provided in the syntenyPlotteR package (https://github.com/Farre-lab/syntenyPlotteR/blob/9623174a881f8e15b19362d4a5985d5e0109b52b/R/DrawSyntenyPlot.R#L212)

draw_linear_laurie <-
  function (output,
             ,
            ...,
            directory = NULL,
            fileformat = "png",
            colours = colours.default,
            w = 13,
            h = 5,
            opacity = 0.5) {
    if (is.null(directory)) {
      directory <- tempdir()
    }
    synteny.data.reframing <-
      function(data, tar.y, ref.y, compiled.size) {
        synteny <- data.frame()
        for (i in c(1:nrow(data))) {
          reference <- data[i, "ref.species"]
          target <- data[i, "tar.species"]
          tar_chr <- data[i, "tarchr"]
          ref_chr <- data[i, "refchr"]
          dir <- data[i, "dir"]
          tar_sizes <- compiled.size[compiled.size$species ==
                                       target,]
          names(tar_sizes) <- c("tarchr", "size", "species",
                                "xstart", "xend")
          ref_sizes <- compiled.size[compiled.size$species ==
                                       reference,]
          names(ref_sizes) <- c("refchr", "size", "species",
                                "xstart", "xend")
          tar_add <- tar_sizes[as.character(tar_sizes$tarchr) ==
                                 as.character(tar_chr),]$xstart
          ref_add <- ref_sizes[as.character(ref_sizes$refchr) ==
                                 as.character(ref_chr),]$xstart
          tar_y <- tar.y
          ref_y <- ref.y
          tar_xstart <- data[i, "tarstart"] + tar_add
          tar_xend <- data[i, "tarend"] + tar_add
          ref_xstart <- data[i, "refstart"] + ref_add
          ref_xend <- data[i, "refend"] + ref_add
          inverted <- grepl("-", dir, fixed = TRUE)
          if (inverted == TRUE) {
            df <- data.frame(
              x = c(tar_xstart, tar_xend,
                    ref_xstart, ref_xend),
              y = c(tar_y, tar_y,
                    ref_y, ref_y),
              fill = ref_chr,
              group = paste0("s",
                             i),
              ref = reference,
              tar = target
            )
          }
          else {
            df <- data.frame(
              x = c(tar_xstart, ref_xstart,
                    ref_xend, tar_xend),
              y = c(tar_y, ref_y, ref_y,
                    tar_y),
              fill = ref_chr,
              group = paste0("s",
                             i),
              ref = reference,
              tar = target
            )
          }
          synteny <- rbind(synteny, df)
        }
        return(synteny)
      }
    colours.default <-
      c(
        `1` = "#BFD73B",
        `2` = "#39ACE2",
        `3` = "#F16E8A",
        `4` = "#2DB995",
        `5` = "#855823",
        `6` = "#A085BD",
        `7` = "#2EB560",
        `8` = "#D79128",
        `9` = "#FDBB63",
        `10` = "#AFDFE5",
        `11` = "#BF1E2D",
        `12` = "purple4",
        `13` = "#B59F31",
        `14` = "#F68B1F",
        `15` = "#EF374B",
        `16` = "#D376FF",
        `17` = "#009445",
        `18` = "#CE4699",
        `19` = "#7C9ACD",
        `20` = "#84C441",
        `21` = "#404F23",
        `22` = "#607F4B",
        `23` = "#EBB4A9",
        `24` = "#F6EB83",
        `25` = "#915F6D",
        `26` = "#602F92",
        `27` = "#81CEC6",
        `28` = "#F8DA04",
        `29` = "peachpuff2",
        `30` = "gray85",
        `33` = "peachpuff3",
        W = "#9590FF",
        Z = "#666666",
        Y = "#9590FF",
        X = "#666666",
        LGE22 = "grey",
        LGE64 = "gray64",
        `1A` = "pink",
        `1B` = "dark blue",
        `4A` = "light green",
        Gap = "white",
        LG2 = "black",
        LG5 = "#CC99CC"
      )
    xstart <-
      xend <-
      refchr <-
      tarchr <- x <- y <- group <- fill <- chromosome <- species <- NULL
    sizes <- utils::read.delim(sizefile, header = FALSE)
    names(sizes) <- c("chromosome", "size", "species")
    sizes$size <- as.numeric(gsub(",", "", sizes$size))
    count <- 0
    compiled.size <- data.frame()
    for (i in unique(sizes$species)) {
      size.intermediate <- sizes[sizes$species == i,]
      for (x in c(1:nrow(size.intermediate))) {
        if (x == 1) {
          total_start <- 1
          total_end <- size.intermediate[x, "size"]
        }
        else {
          total_start <- total_end + 6e+06
          total_end <- total_start + size.intermediate[x,
                                                       "size"]
        }
        size.intermediate[x, "xstart"] <- total_start
        size.intermediate[x, "xend"] <- total_end
      }
      compiled.size <- rbind(compiled.size, size.intermediate)
    }
    for (z in unique(compiled.size$species)) {
      compiled.size$y[compiled.size$species == z] <- count
      count <- count + 2
    }
    list.of.files <- list()
    for (i in list(...)) {
      list.of.files[[i]] <- i
    }
    listsynt <- list()
    for (i in 1:length(list.of.files)) {
      num <- i
      file <- list.of.files[[num]]
      dataTMP <- utils::read.delim(file, header = FALSE)
      data2 <- dataTMP[, c(4, 5, 6, 1, 2, 3, 7, 8, 9)]
      colnames(data2) <- c(
        "tarchr",
        "tarstart",
        "tarend",
        "refchr",
        "refstart",
        "refend",
        "dir",
        "ref.species",
        "tar.species"
      )
      data2$tarstart <- as.numeric(gsub(",", "", data2$tarstart))
      data2$tarend <- as.numeric(gsub(",", "", data2$tarend))
      data2$refstart <- as.numeric(gsub(",", "", data2$refstart))
      data2$refend <- as.numeric(gsub(",", "", data2$refend))
      reference <- data2[1, "ref.species"]
      target <- data2[1, "tar.species"]
      ref_y <- compiled.size[compiled.size$species == reference,
                             "y"]
      tar_y <- compiled.size[compiled.size$species == target,
                             "y"]
      if (tar_y[1] > ref_y[1]) {
        ref_y <- ref_y[1] + 0.1
        tar_y <- tar_y[1]
      }
      else {
        ref_y <- ref_y[1]
        tar_y <- tar_y[1] + 0.1
      }
      x <- synteny.data.reframing(data2, tar_y, ref_y, compiled.size)
      x$fill <- as.factor(x$fill)
      listsynt[[i]] <- x
    }
    compiled.size$chromosome <- as.factor(compiled.size$chromosome)
    
    p <- ggplot2::ggplot()
    for (i in 1:length(listsynt)) {
      data <- listsynt[[i]]
      reference <- data[1, "ref"]
      target <- data[1, "tar"]
      ref_sizes <- compiled.size[compiled.size$species == reference,]
      tar_sizes <- compiled.size[compiled.size$species == target,]
      p <-
        p +
        ggplot2::geom_rect(
          data = ref_sizes,
          mapping = ggplot2::aes(
            xmin = xstart,
            xmax = xend,
            ymin = y,
            ymax = y + 0.1,
            fill = chromosome
          ),
          color = "black",
          alpha = 0.85,
          linewidth = 0.1
        ) +
        
        ggplot2::geom_polygon(
          data = data,
          alpha = opacity,
          ggplot2::aes(
            x = x,
            y = y,
            group = group,
            fill = fill
          )
        ) +
        # Ref chr tags
        ggplot2::geom_text(
          data = ref_sizes,
          ggplot2::aes(
            x = (xstart + xend) / 2,
            y = y+0.2,
            label = chromosome
          ),
          size = 3,
          angle = 0,
        ) +
        
        # Ref pecies tag
        ggplot2::geom_text(
          data = ref_sizes,
          mapping = ggplot2::aes(x = 2, y = y, label = species),
          size = 3,
          hjust = 1,
          vjust = 1
        ) +
        # Query chr rectangles
        ggplot2::geom_rect(
          data = tar_sizes,
          mapping = ggplot2::aes(
            xmin = xstart,
            xmax = xend,
            ymin = y,
            ymax = y + 0.1
          ),
          fill = "grey90",
          color = "black",
          alpha = 0.85,
          linewidth = 0.1
        ) +
        # Query chr rectangles
        ggplot2::geom_text(
          data = tar_sizes,
          ggplot2::aes(
            x = (xstart + xend) / 2,
            y = y + 0.05,
            label = chromosome
          ),
          size = 2,
          angle = 45
        ) +
        # Query species tag
        ggplot2::geom_text(
          data = tar_sizes,
          mapping = ggplot2::aes(x = 2, y = y, label = species),
          size = 3,
          hjust = 1)
          
        
    }
    p <- p +
      ggplot2::scale_fill_manual(values = colours) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      )
    message(paste0("Saving linear image to ", directory))
    print(p)
    ggplot2::ggsave(
      paste0(directory, "/", output, ".", fileformat),
      p,
      device = fileformat,
      width = w,
      height = h
    )
  }
