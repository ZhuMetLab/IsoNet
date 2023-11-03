# Copyright © 2022- ZhuMSLab ALL right reserved
#
# @author: Yan-dong Yin
# @contact: yddream@gmail.com
# @project: MIDNet
# @file: zzz.R
# @time: 2023-04-14 22:51


################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  nmark <- 80
  welcome <- paste0("WELCOME to ", getPackageName(), " v", packageVersion(getPackageName()))
  welcome <- .gen_marks(welcome)
  end_mark <- .gen_marks(nmark = nchar(welcome))
  packageStartupMessage(welcome,
                        "\nIf you have any questions, please send email to zhouzw@sioc.ac.cn.",
                        "\nAuthors: Yang Gao, Mingdu Luo, Hongmiao Wang, Zhiwei Zhou and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).",
                        "\nMaintainer: Yang Gao",
                        "\n",.gen_marks(content = 'NEWS', mark='-', nmark= nmark),
                        "\nVersion 0.2.13 (20231020)",
                        "\n\to Debug removing features without MS2 during file checking.",
                        "\n",
                        end_mark,
                        "\n")
}

.gen_marks <- function(content=NULL, mark='=', nmark= 80) {
  if (!is.null(content)) {
    mark_content <-  paste0(rep(mark, (nmark - nchar(content))/2-1), collapse = '')
    return(paste(mark_content, content, mark_content))
  } else {
    return(paste0(rep(mark, nmark), collapse = ''))
  }
}