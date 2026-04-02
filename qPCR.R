library(shiny)
library(DT)
library(tidyverse)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(car)
library(PMCMRplus)
library(openxlsx)
library(colourpicker)
library(ggsci)
library(ggpubr)

#函数#####
##计算最终CT值得函数#####
calculate_expression <- function(df, 
                                 group,
                                 reference_group,
                                 sample,
                                 gene,
                                 reference_gene ,
                                 CT
                                 ) {
  
  df[[CT]]<- as.numeric(df[[CT]])
  # 1. 计算 CT 均值 (技术重复)
  df_summarized=df|>
    group_by(across(all_of(c(gene,sample,group))))|>
    summarise(CT_mean = mean(.data[[CT]], na.rm = TRUE), .groups = "drop")|>
    drop_na()
  
  # 计算内参CT值均数
  A <- df_summarized |>
    filter(.data[[gene]] == reference_gene) |>
    group_by(across(all_of(c(group, sample)))) |>
    summarise(CT_mean_reference_gene = round(mean(CT_mean), 6),.groups = 'drop')
  
  df_delta <- df_summarized |>
    left_join(A, by =c( sample,group))
  
  # 计算ΔCT = 目的基因CT - 内参CT均值
  df_delta <- df_delta |>
    mutate(deltaCT = CT_mean - CT_mean_reference_gene)
  # 计算目的基因对照组ΔCT均值
  B <- df_delta |>
    filter(.data[[gene]] != reference_gene, .data[[group]] == reference_group) |>
    group_by(.data[[gene]]) |>
    summarise(mean = mean(deltaCT))
  
  df_delta <- df_delta |>
    left_join(B, by = gene )|>
    drop_na()

    
  
  # 计算ΔΔCT = 目的基因ΔCT - 对照组ΔCT
  df_delta <- df_delta |>
    mutate(deltadeltaCT = deltaCT - mean)
  
  # 计算2^-ΔΔCT
  final_df <- df_delta |>
    mutate(result = 2^(-deltadeltaCT))|>
    arrange(.data[[gene]],.data[[group]])
  
  return(final_df)
}

##正态性检验#####
zhengtai=function(df,group,test_gene,gene){
  if(test_gene == 'ALL'){
  df <- df |> mutate(result = as.numeric(result))
  df|>
    group_by(across(all_of(c(group,gene))))|>
    summarise(
      statistics=shapiro.test(result)$statistic,
      PValue=shapiro.test(result)$p.value)|>
      mutate(is_normal = ifelse(PValue > 0.05, "是", "否"))
  }else{
    df|>
      filter(.data[[gene]] == test_gene)|>
      group_by(across(all_of(group)))|>
      summarise(
        statistics=shapiro.test(result)$statistic,
        PValue=shapiro.test(result)$p.value)|>
      mutate(is_normal = ifelse(PValue > 0.05, "是", "否"))
    
    }
}

##方差齐性#####
levene_test <- function(df,group,test_gene,gene_col){
  sub_df <- df[df[[gene_col]] == test_gene, ]
  
  if (nrow(sub_df) < 2 || length(unique(sub_df[[group]])) < 2) {
    return('数据量过少，暂不支持分析') # 或者返回一个自定义的错误信息，避免后续计算崩溃
  }
  
  res <- tryCatch({
    leveneTest(
      y=sub_df[['result']],
      group = as.factor(sub_df[[group]])
    )
  }, error = function(e) {
    return(NULL)
  })
  
  return(res)
}

##事后检验转换函数#####
format_posthoc_res <- function(res_obj) {
  # 1. 提取 P 值矩阵
  p_mat <- res_obj$p.value
  
  # 2. 将下三角矩阵转为长表格
  # 使用 as.table 会自动忽略 NA (上三角部分)
  df <- as.data.frame(as.table(p_mat))
  colnames(df) <- c("组别_1", "组别_2", "P_Value")
  
  # 3. 过滤掉自身对比或重复的行
  df <- df[!is.na(df$P_Value), ]
  
  # 4. 增加显著性星号标记
  df$显著性 <- cut(df$P_Value, 
                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                labels = c("***", "**", "*", "ns"))
  
  # 5. 格式化 P 值显示 (保留 4 位小数)
  df$P_Value <- format.pval(df$P_Value, digits = 4, eps = 0.001)
  
  return(df)
}

##组间分析#####
between_group_test <- function(df,method,test_gene,gene_col,group){
  #筛选出特定基因的数据
  sub_df <- df[df[[gene_col]] == test_gene, ]
  
  if(method == 'ANOVA'){
    formula_str <- as.formula(sprintf("result ~ `%s`", group))
    fit <- aov(formula_str,data = sub_df)
    return(summary(fit))
  }else if(method == 'Kruskal-Wallis test'){
    formula_str <- as.formula(sprintf("result ~ `%s`", group))
    fit <- kruskal.test(formula_str, data = sub_df)
    return(fit)
  }

}

##绘图
draw_plot <- function(data,
                      x_var,
                      y_var, 
                      plot_type = 'box',
                      palette = 'npg',
                      bin_width = 0.8,                    
                      xlabs = 'Group',
                      ylabs = 'Expression',                      
                      font_size = 14,
                      font_family = 'serif',
                      alpha = 1,
                      show_points = T,
                      show_points_color = 'black'
){

  # 基础图形构建
  p <- data|>
    ggplot(aes(x=.data[[x_var]],
               y=.data[[y_var]],
               fill = .data[[x_var]]))
  if (plot_type == "box"){
    p <- p + geom_boxplot(width = bin_width, alpha = alpha, outlier.shape = NA)
  }else{
    p <- p + stat_summary(fun = mean, geom = "bar", width = bin_width, alpha = alpha, color = "black") +
      stat_summary(fun.data = mean_se, fun.args = list(mult = 1), 
                   geom = "errorbar", width = 0.2)
  }
  # 抖动点 (Jitter points)
  if (show_points) {
    p <- p + geom_jitter(width = 0.15, alpha = 0.5, shape = 21, fill = show_points_color)
  }
  
  ## 配色方案 (ggsci)
  p <- p + switch(palette,
                  "npg"    = ggsci::scale_fill_npg(),
                  "aaas"   = ggsci::scale_fill_aaas(),
                  "nejm"   = ggsci::scale_fill_nejm(),
                  "lancet" = ggsci::scale_fill_lancet(),
                  "jco"    = ggsci::scale_fill_jco(),
                  'ucscgb' = ggsci::scale_fill_ucscgb(),
                  'd3'     = ggsci::scale_fill_d3(),
                  'flatui' = ggsci::scale_fill_flatui(),
                  ggsci::scale_fill_npg()
  )
  
  ## 主题与字体
  p <- p+theme_classic()+
    labs(x=xlabs,y=ylabs)+
    theme(
      text = element_text(family = font_family,size = font_size),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      legend.position = 'none'
    )+scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  return(p)
}



#UI#####
##第一页UI#####
calc_ui <- if(T){fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      shinyFeedback::useShinyFeedback(),
      width = 3,
      fileInput('file','上传csv文件',accept = '.csv',width = '100%'),
      
      selectInput('group', '请选择实验分组', choices = NULL),
      
      conditionalPanel(
        condition = "input.group != ''",
        selectInput('reference_group', '选择对照组', choices = NULL)
      ),
      
      selectInput('gene', '请选择基因列', choices = NULL),
      
      conditionalPanel(
        condition = "input.gene != ''",
        selectInput('reference_gene', '选择内参基因', choices = NULL)
      ),
      
      selectInput('sample', '请选择样本列', choices = NULL),
      selectInput('CT', '请选择CT值列', choices = NULL),

      hr(),
      splitLayout(
        cellWidths = c("50%", "50%"),
        actionButton(
          'Start_calc', 
          '开始计算', 
          class = "btn-primary", 
          style = "width: 100%"
        ),
        downloadButton(
          'download', 
          '导出结果', 
          class = "btn-success", 
          style = "width: 100%"
        )
      ),
      br(),
      splitLayout(
        cellWidths = c('50%','50%'),
        actionButton('Run_demo','运行内置数据',
                     icon = icon('play-circle'),
                     class = 'btn-info',
                     style = 'width:100%'),
        downloadButton('download_manual',
                       '使用说明',
                       class = "btn-default",
                       style = "width: 100%"
        )
      )
    ),
    
    mainPanel(
      width = 9,
      uiOutput('result_box'),

      conditionalPanel(
        condition = "input.Start_calc > 0", # 只有点击开始计算后才显示
        fluidRow(
          box(
            title = tagList(icon("file-alt"), " 分析报告摘要"),
            status = "info",
            solidHeader = TRUE,
            width = 12,
            htmlOutput("analysis_report")
          )
        )
      ),
    )
  )
)
}

##正态性检验UI#####
norm_ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput('normalize_test_gene',label = '请选择正态性检验的基因',choices = c('所有基因')),
      splitLayout(
        cellWidths = c("50%", "50%"),
        actionButton(
          'start_normalize_test', 
          '开始分析', 
          class = "btn-primary", 
          style = "width: 100%"
        ),
        downloadButton(
          outputId = 'download_normalize_test_result', 
          label = '导出结果', 
          class = "btn-success", 
          style = "width: 100%"
        )
      ),
    ),
    mainPanel(
      width = 9,
      box(title = '正态分析结果',
          width = 12,
          solidHeader = T,
          status = 'primary',
          DTOutput('normalize_test_result')
      ),
    )
  )
)


##方差齐性检验UI#####
levene_ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput('levene_test_gene',label = '请选择方差齐性检验的基因',choices = ''),
      splitLayout(
        cellWidths = c("50%", "50%"),
        actionButton(
          'start_levene_test', 
          '开始分析', 
          class = "btn-primary", 
          style = "width: 100%"
        ),
      ),
    ),
    mainPanel(
      width = 9,
      box(
        title = '方差齐性检验结果',
        width = 12,
        solidHeader = T,
        status = 'primary',
        verbatimTextOutput('levene_test_result')
      )
    )
  )
)

##组间分析UI#####
between_group_test_ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput('between_group_test_method',
                  label = '请选择组间检验方法',
                  choices = c('方差分析（ANOVA）'='ANOVA','秩和检验（KW H检验）'='Kruskal-Wallis test')),
      selectInput('between_group_test_gene',
                  label= '请选择组间检验的基因',
                  choices = ''),
      splitLayout(
        cellWidths = c('50%','50%'),
        actionButton('start_between_group_test',
                     label ='开始分析', 
                     class = "btn-primary",
                     style = "width: 100%")
      )
      
    ),
    mainPanel(
      width = 9,
      conditionalPanel(
        condition = 'input.start_between_group_test > 0',
        fluidRow(
          box(
            title = '组间分析结果',
            status = 'primary',
            solidHeader = T,
            width = 12,
            verbatimTextOutput('between_group_test_result')
          )
        )
      )
    )
  )
)

##事后检验UI#####
post_hoc_ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput('post_hoc_gene','请选择事后检验的基因',choices = ''),
      pickerInput(
        inputId = "post_hoc_methods",
        label = "选择事后检验方法:",
        choices = list(
          "参数检验 (ANOVA 后)" = c("Tukey HSD" = "Tukey", 
                                    "Dunnett" = "Dunnett",
                                    'Bonferroni' = 'Bonferroni',
                                    "LSD" = "LSD"),
          "非参数检验 (KW H 后)" = c("Dunn's Test" = "Dunn",
                               "Wilcoxon Pairwise" = "Wilcoxon",
                               'Nemenyi Test'='Nemenyi')),                               
        selected = "Tukey",
        options = list(
          `actions-box` = TRUE,      # 开启全选/取消全选按钮
          `deselect-all-text` = "全不选", # 自定义按钮文字
          `select-all-text` = "全选",     # 自定义按钮文字
          `none-selected-text` = "未选择任何方法"),# 未选中时的提示
        multiple = TRUE
      ),
      splitLayout(
        cellWidths = c('50%','50%'),
        actionButton(
          inputId = 'start_post_hoc',
          label = '开始事后检验',
          class = "btn-primary",
          style = "width: 100%;"),
        
        downloadButton(
          outputId = "download_all_posthoc",
          label = "导出所有结果",
          class = "btn-success", 
          style = "width: 100%;")
      ),
      
    ),
    mainPanel(
      width = 9,
      box(
        title = '各种事后检验方法说明',
        status = 'primary',
        solidHeader = T,
        width = 12,
        uiOutput('Explanation_text'),
        collapsible = T,
      ),
      # 动态生成的结果区域
      uiOutput('dynamic_post_hoc_results')
      
    )
  )
)

##绘图ui#####
plot_ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h4("1. 基础设置", style = "font-weight: bold; color: #3c8dbc;"),
      selectInput("plot_gene", "选择分析基因:", choices = ""),
      
      selectInput("plot_type", "图形种类:", 
                   choices = c("箱线图" = "box", "柱状图" = "bar"),
                   ),
      uiOutput("group_order_ui"),
      
      hr(),
      tags$h4("2. 统计标注", style = "font-weight: bold; color: #3c8dbc;"),
      checkboxInput("show_points", "显示原始数据点", value = FALSE),
      checkboxInput("show_p_val", "显示显著性标记 (该功能暂未实现)", value = FALSE),
      
      
      hr(),
      actionButton("render_plot", "生成图片", class = "btn-success", width = "100%")
    ),
    mainPanel(
      width = 9,
      # 图片预览区
      box(
        title = '图片预览',
        status = 'primary',
        solidHeader = T,width = 12,
        plotOutput('main_plot',height = '500px')
      ),
      #自定义参数
      tabBox(
        title = tagList(icon('gear'),'自定义细节'),
        id = "plot_tabs", width = 12,
        tabPanel("颜色与样式",
                 fluidRow(
                   column(4,
                          pickerInput(inputId = 'plot_palette',
                                      label = '选择分组配色',
                                      choices = list(
                                          "常用经典" = c("NPG (Nature)" = "npg", 
                                                     "AAAS (Science)" = "aaas", 
                                                     "NEJM (医学经典)" = "nejm"),
                                          "渐变/鲜艳" = c("Lancet (柳叶刀)" = "lancet", 
                                                      "JCO (临床肿瘤)" = "jco", 
                                                      "UCSC (染色体色)" = "ucscgb"),
                                          "扁平化" = c("D3 (数据可视化)" = "d3", 
                                                    "Flat UI" = "flatui")
                                        ),
                                      selected = "npg",
                                      options = list(`live-search` = TRUE)
                                      ),
                          ),
                   column(4,sliderInput('bin_width',label = '柱子/箱体宽度',min = 0.2,max = 1,value = 0.8)),
                   column(4,conditionalPanel(
                     condition = 'input.show_points == true',
                     pickerInput(inputId = 'show_points_color',
                                 label = '数据原点颜色',
                                 choices = c('black','white','yellow','red','blue'),selected = 'black')
                   )),
                   
                 ),
                 fluidRow(
                   column(12, 
                          checkboxInput("fill_transparent", "半透明填充 ", value = FALSE),
                          conditionalPanel(
                            condition = "input.fill_transparent == true",
                            sliderInput("plot_alpha", "颜色透明度:", min = 0.1, max = 1, value = 0.7)
                          )
                   )
                 ),
                 ),
        tabPanel('坐标轴与文字',
                 fluidRow(
                   column(4,textInput('xlab','X轴标签',value = 'Group')),
                   column(4,textInput('ylab','y轴标签',value = 'Expression')),
                   column(4, sliderInput("font_size", "全局字体大小:", min = 8, max = 24, value = 14)),
                   column(4,selectInput('font_family','字体样式',
                                        choices = c('Times New Roman (有衬线)' = 'serif',
                                                    'Arial (无衬线)' = 'sans',
                                                    'Courier New (等宽)' = 'mono'),
                                        selected = 'serif'))
                 )),
        tabPanel("导出设置",
                 fluidRow(
                   column(4, selectInput("file_type", "文件格式:", choices = c("pdf", "png", "tiff"))),
                   column(4, numericInput("plot_res", "分辨率 (DPI):", value = 300)),
                   column(4, numericInput('plot_width','图片宽度(英寸)',value = 8)),
                   column(4, numericInput('plot_height','图片高度(英寸)',value = 6)),
                   column(4, style = "margin-top: 25px;", downloadButton("download_plot", "下载图片"))
                   )),
      ),
    )
  )
) 




##联系UI#####
contact_ui <- fluidPage(
  br(),
  column(width = 8, offset = 2,
         # 1. 工具简介
         wellPanel(
           style = "background: #fff; border-top: 3px solid #3c8dbc;",
           tags$h2("qPCR 数据自动化分析系统", style = "color: #3c8dbc; font-weight: bold;"),
           tags$p("本工具旨在为生物医学研究者提供一站式的 qPCR 原始数据处理与统计分析。", 
                  style = "font-size: 16px; color: #666;"),
           hr(),
           
           # 2. 核心功能点
           tags$h4("主要功能：", style = "font-weight: bold;"),
           tags$ul(
             tags$li("自动化 2^-ΔΔCt 计算流程"),
             tags$li("支持统计分析及多种事后检验"),
             tags$li("内置 Nature/Science/Lancet 等顶级期刊配色"),
             tags$li("导出高分辨率 (300+ DPI) 的 PDF/PNG 矢量图")
           )
         ),
         
         # 3. 作者与技术支持（第一行信息卡片）
         fluidRow(
           column(6,
                  box(
                    title = "技术支持", status = "primary", solidHeader = TRUE, width = 12,
                    tags$p(icon("envelope"), " 联系邮箱：3170070056@smu.deu.cn"),
                    tags$p(icon("github"), " GitHub: ", tags$a(href="https://github.com/Liqiuhuang", "https://github.com/Liqiuhuang", target="_blank"))
                  )
           ),
           column(6,
                  box(
                    title = "版本信息", status = "info", solidHeader = TRUE, width = 12,
                    tags$p("当前版本：v2.0.1"),
                    tags$p("更新日期：2026-03"),
                   
                  )
           )
         ),
         
         # 4. 二维码展示区（第二行图片对齐）
         fluidRow(
           style = "margin-top: 20px;",
           column(6,
                  div(style = "text-align: center;",
                      img(src = "2.jpg", 
                          height = "250px", 
                          style = "border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);"),
                      p(style = "margin-top: 15px; color: #666; font-style: italic;", "如有问题，联系作者")
                  )
           ),
           column(6,
                  div(style = "text-align: center;",
                      img(src = "1.jpg", 
                          height = "250px", 
                          style = "border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);"),
                      p(style = "margin-top: 15px; color: #666; font-style: italic;", "制作不易，支持作者")
                  )
           )
         )
  )
)


##总UI#####
ui <- dashboardPage(
  dashboardHeader(title='qPCR 数据分析系统'),
  dashboardSidebar(
    sidebarMenu(
      menuItem('数值计算',tabName = 'calc',icon = icon('calculator')),
      menuItem('统计分析',tabName = 'stat',icon = icon('th'),
               menuSubItem('正态检验',tabName = 'normality'),
               menuSubItem('方差齐性检验',tabName = 'levene'),
               menuSubItem('组间分析',tabName = 'between_group_test'),
               menuSubItem('事后检验',tabName = 'post-hoc')),
      menuItem('图片绘制',tabName = 'figure',icon = icon('chart-bar')),
      menuItem('联系作者',tabName = 'contact',icon = icon('envelope'))
    )
  ),
  dashboardBody(
    if(T){
    tags$head(
      tags$style(HTML("
      /* 1. 锁定 Header (顶部栏) */
      .main-header {
        position: fixed;
        width: 100%;
        top: 0;
        left: 0;
        z-index: 1030; /* 确保在最上层 */
      }
      
      /* 2. 锁定 Sidebar (侧边栏) */
      .main-sidebar {
        position: fixed !important;
        top: 0;
        left: 0;
        padding-top: 50px !important; /* 预留出 Header 的高度 */
        height: 100vh !important;
        overflow-y: auto; /* 如果侧边栏菜单太长，允许内部滚动 */
      }
      
      /* 3. 修正主体内容区域的偏移 */
      .content-wrapper {
        padding-top: 50px !important; /* 确保不被固定的 Header 遮挡 */
        min-height: 100vh !important;
      }
      
      /* 4. 解决侧边栏折叠后的宽度问题 (可选优化) */
      .sidebar-mini.sidebar-collapse .main-header .logo { width: 50px !important; }
      .sidebar-mini.sidebar-collapse .main-header .navbar { margin-left: 50px !important; }
      .sidebar-mini.sidebar-collapse .content-wrapper { margin-left: 50px !important; }
    "))
    )},
    
    tabItems(
      #第一页布局
      tabItem(tabName = 'calc',calc_ui),
      #第二页布局
      tabItem(tabName = 'stat'),
      tabItem(tabName = 'normality',norm_ui),
      tabItem(tabName = 'levene',levene_ui ),
      tabItem(tabName = 'between_group_test',between_group_test_ui),
      tabItem(tabName = 'post-hoc',post_hoc_ui),
      #第三页布局
      tabItem(tabName = 'figure',plot_ui),
      #第四页
      tabItem(tabName = 'contact',contact_ui)
    )
    
  )
)

#server#####

server <- function(input, output, session) {
  
  ##第一页server#####
  
  df <- reactiveVal(NULL)
  #读取数据
  observeEvent(input$file,{    
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    is_csv <- ext == 'csv'
    shinyFeedback::feedbackDanger(
      'file',
      !is_csv,
      '格式错误：请上传csv格式文件'
    )
    validate(
      need(is_csv, "请务必上传 csv 格式文件")
    )
    data <- vroom::vroom(file = input$file$datapath,delim = ',')
    df(data)
  })
  
  #当文件上传后，更新所有列名相关的下拉框
  observeEvent(df(), {
    cols <- c("请选择..." = "", colnames(df()))
    updateSelectInput(session, "group", choices = cols)
    updateSelectInput(session, "sample", choices = cols)
    updateSelectInput(session, "gene", choices = cols)
    updateSelectInput(session, "CT", choices = cols)
  })
  
  #当分组列确定后，更新对照组选项
  observeEvent(input$group, {
    req(input$group != "")
    choices <- unique(df()[[input$group]])
    updateSelectInput(session, "reference_group", choices = c("请选择..." = "", choices))
})
  
  #当基因列确定后，更新内参基因选项
  observeEvent(input$gene, {
    req(input$gene != "")
    choices <- unique(df()[[input$gene]])
    updateSelectInput(session, "reference_gene", choices = c("请选择..." = "", choices))
  })
  
  
  #计算最终结果
  result <- eventReactive(input$Start_calc, {
    # 门禁：确保这些输入都不为空
    req(
        input$group, 
        input$reference_group, 
        input$sample, 
        input$gene, 
        input$reference_gene, 
        input$CT)
    
    calculate_expression(
      df = df(),
      group = input$group,
      reference_group = input$reference_group,
      sample = input$sample,
      gene = input$gene,
      reference_gene = input$reference_gene,
      CT = input$CT
    )
    
    })
  
  #动态渲染结果Box
  output$result_box <- renderUI({
    req(result())
    fluidRow(
      box(
        title = '分析后的数据',
        status = 'primary',
        solidHeader = T,
        DTOutput('data_table'),
        width = 12
      )
    )
  })
  
  #渲染结果box里面的DT
  output$data_table <- renderDT({
    req(result())
    
    datatable(
      result(),
      filter = 'top',
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Chinese.json')
      )
    ) %>% 
      # 将数值列格式化为 6 位小数
      formatRound(columns = c('CT_mean','CT_mean_reference_gene','deltaCT','mean','deltadeltaCT','result'), digits = 6)
  })
  
  #渲染分析报告摘要
  output$analysis_report <- renderUI({
    # 确保已经有了计算结果
    req(result())
    
    # 提取一些统计信息
    gene_list <- unique(result()[[input$gene]])
    sample_count <- length(unique(result()[[input$sample]]))
    group_list <- unique(result()[[input$group]])
    
    # 将列表转换为带样式的 HTML 字符串
    gene_tags <- paste0("<span class='label label-default' style='margin: 0 3px; display:inline-block; font-size: 14px;'>", 
                        gene_list, "</span>", collapse = " ")
    group_tags <- paste0("<span class='label label-default' style='margin: 0 3px; display:inline-block; font-size: 14px;'>", 
                         group_list, "</span>", collapse = " ")
    
    # 使用 HTML 拼接字符串
      HTML(paste0(
        "<div style='line-height: 2.5; font-size: 16px; padding: 20px; color: #333;'>",

        
        # 第一行：基因信息
        "<b>● 基因情况：</b>本次分析共包含 <b>", length(gene_list), "</b> 个基因，分别为：", 
        "<span>", gene_tags, "</span><br/>",
        
        # 第二行：分组信息
        "<b>● 实验分组：</b>划分为 <b>", length(group_list), "</b> 个实验分组，分别为：", 
        "<span>", group_tags, "</span><br/>",
        
        # 第三行：样本信息
        "<b>● 样本规模：</b>涉及 <b>", sample_count, "</b> 个样本（已完成技术重复均值计算）。<br/>",
        
        "<h3 style='color: #2c3e50; border-bottom: 2px solid #3c8dbc; padding-bottom: 8px; margin-top: 25px;'>结果列名说明</h3>",
        "<div style='background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 5px solid #3c8dbc;'>",
        "<ul style='list-style-type: none; padding-left: 5px; line-height: 2;'>",
        "<li> <b>CT_mean</b>：技术性重复孔的 CT 原始值均值。</li>",
        "<li> <b>CT_mean_reference_gene</b>：该样本内参基因（", input$reference_gene, "）的 CT 均值。</li>",
        "<li> <b>ΔCT</b> = CT_mean - CT_mean_reference_gene</li>",
        "<li> <b>mean</b>：目的基因对照组 ΔCT 的平均值。</li>",
        "<li> <b>ΔΔCT</b> = ΔCT - mean</li>",
        "<li> <b>最终结果 (result)</b> = 2<sup>-ΔΔCT</sup></li>",
        "</ul>",
        "</div>",
    
        "<h4 style='color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 5px; margin-top: 15px;'>结果提示</h4>",
        "<span style='color: #e67e22;'>⚠️ 注意：已从结果列表中剔除了内参基因数据。</span>",

        "</div>"
    ))
  })
  
  #渲染下载按钮
  output$download <- downloadHandler(
    filename = function() {
      # 去掉原始扩展名再拼接，避免出现 .csv.csv
      original_name <- tools::file_path_sans_ext(input$file$name)
      paste0('分析结果_', original_name, '_', Sys.Date(), '.csv')
    }, 
    content = function(file) {
      # row.names = FALSE 可以防止生成的csv多出一列无用的行号
      write.csv(result(), file, row.names = FALSE)
    }
  )
  
  #加载内置数据
  observeEvent(input$Run_demo,{
    if(file.exists('test.csv')){
      demo_data <- read.csv('test.csv',stringsAsFactors = FALSE)
      
      df(demo_data)
      shinyFeedback::showToast("success", "已成功加载内置测试数据！")
    }
  })
  #下载使用说明
  output$download_manual <- downloadHandler(
  filename = function() { "使用说明.docx" },
  content = function(file) {
    # 确保 www 文件夹下有这个文件，或者直接指定路径
    file.copy("使用说明.docx", file)
  }
  )
  
  
  
  
  ##正态性检验server#####
  #更新正态性检验的基因选项
  observeEvent(
    req(result()),
    updateSelectInput(inputId = 'normalize_test_gene',choices = c('所有基因'='ALL',unique(result()[[input$gene]])))
  )
  
  #进行正态性检验
  normalize_test_result_data <- eventReactive(input$start_normalize_test,{
    req(result(),input$group,input$gene)
    zhengtai(df = result(),
             group = input$group,
             test_gene = input$normalize_test_gene,
             gene = input$gene)

  })
  
  #渲染正态性结果    
  output$normalize_test_result <-renderDT({
      req(normalize_test_result_data())
      
      datatable(
        normalize_test_result_data(),
        filter = 'top',
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          autoWidth = FALSE,
          columnDefs = list(list(className = 'dt-center', targets = "_all")),
          language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Chinese.json')
        )
      ) %>% 
        # 将数值列格式化为 6 位小数
        formatRound(columns = c('statistics','PValue'), digits = 6)
    })
  
  #渲染下载正态性检验结果
  output$download_normalize_test_result <- downloadHandler(
    filename = function() {
      # 去掉原始扩展名再拼接，避免出现 .csv.csv
      original_name <- tools::file_path_sans_ext(input$file$name)
      paste0('正态性检验分析结果_', original_name, '_', Sys.Date(), '.xlsx')
    }, 
    content = function(file) {
      writexl::write_xlsx(normalize_test_result_data(), path = file)
    }
  )
  
  ##方差齐性server#####
  #更新方差齐性检验的基因选项
  observeEvent(result(),{
    req(result(),input$gene)
    updateSelectInput(inputId = 'levene_test_gene',choices = unique(result()[[input$gene]]))})
  
  #进行方差齐性检验
  levene_test_result <- eventReactive(input$start_levene_test,{
    req(result(),
        input$group,
        input$gene,
        input$levene_test_gene)
    
    res_obj <- levene_test(df = result(),
                           group = input$group,
                           test_gene = input$levene_test_gene,
                           gene_col = input$gene)
    
    return(list(
      analysis_res = res_obj,
      gene_name = input$levene_test_gene
    ))
    
  })
  
  #渲染方差齐性检验结果
  output$levene_test_result <- renderPrint({
    req(levene_test_result())
    
    full_data <- levene_test_result()
    res <- full_data$analysis_res
    current_gene <- full_data$gene_name
    
    ## 如果 levene_test 函数内部 tryCatch 返回了 NULL，这里要处理
    if (is.null(res)) {
      return(cat("错误：当前数据无法执行 Levene 检验（可能样本量过少或存在大量缺失值）。"))
    }
    
    #提取第一行数据
    p_val  <- res[["Pr(>F)"]][1]
    df_val <- res[["Df"]][1]
    f_val  <- res[["F value"]][1]
    
    # 3. 打印详细信息
    cat("--- 详细统计数据表 ---\n")
    print(res)      
    cat("--------------------------------------------\n")
    cat("--- 方差齐性检验结果 (Levene's Test) ---\n")
    cat("目标基因：", current_gene, "\n")
    cat("自由度 (Df)：", df_val, "\n")
    cat("F 统计量：", round(f_val, 3), "\n")
    cat("P 值：", format.pval(p_val, digits = 4, eps = 0.001), "\n")
    
    cat("--------------------------------------------\n")
    
    if (p_val < 0.05) {
      cat("【结论】 P < 0.05，拒绝原假设。\n")
      cat(sprintf("=> 基因 [%s] 在各组间的方差存在显著差异（方差不齐）。\n", current_gene))
      cat("建议：后续分析请使用非参数检验（Kruskal-Wallis）。")
    } else {
      cat("【结论】 P >= 0.05，接受原假设。\n")
      cat(sprintf("=> 基因 [%s] 在各组间的方差齐同（满足方差齐性）。\n", current_gene))
      cat("建议：可以继续进行标准方差分析 (ANOVA)。")
    }
    
    
    }
  )
  
  
  
  ##组间分析server#####
  #更新组间分析的基因选项
  observeEvent(result(),{
    req(result(),input$gene)
    updateSelectInput(inputId = 'between_group_test_gene',choices = unique(result()[[input$gene]]))}
  )
  
  #进行组间分析
  between_group_test_result <- eventReactive(input$start_between_group_test,{
    req(result(),
        input$group,
        input$gene,
        input$between_group_test_method,
        input$between_group_test_gene)
    res_obj <- between_group_test(df = result(),
                                  method =input$between_group_test_method,
                                  test_gene = input$between_group_test_gene,
                                  gene_col = input$gene,
                                  group = input$group 
                                  )
    return(list(
      analysis_res = res_obj,
      gene_name = input$between_group_test_gene
    ))
    
  })
  
  #渲染组间分析
  output$between_group_test_result <- renderPrint({
    req(between_group_test_result())  
    
    full_data <- between_group_test_result()
    res <- full_data$analysis_res
    current_gene <- full_data$gene_name
    
    if((inherits(res, "summary.aov") || (is.list(res) && !inherits(res, "htest")))){

      #提取第一行数据
      stats_row <- res[[1]]
      p_val <- stats_row[["Pr(>F)"]][1]
      df_val <- stats_row[["Df"]][1]
      f_val <- stats_row[["F value"]][1]
      
      cat("\n详细数据表：\n")
      print(res)      
      cat("------------------------\n")
      cat("--- 组间分析统计结果 ---\n")
      cat("分析方法：方差分析",  "\n")
      cat("组间自由度 (Df)：", df_val, "\n")
      cat("F 统计量：", round(f_val, 3), "\n")
      cat("P 值：", format.pval(p_val, digits = 4, eps = 0.001), "\n")
      if(p_val < 0.05){
        cat(sprintf('结论：可以认为基因 [%s] 在不同分组间的表达水平显著不同。',
                    current_gene))
      }else{
        cat(sprintf('结论：目前没有足够证据表明基因 [%s] 在各组间的表达量存在显著差异。',
                    current_gene))
        
      }
    }else if(inherits(res, "htest")){
      
      p_val  <- res$p.value
      df_val <- res$parameter
      chi_sq <- res$statistic
      
      cat("\n详细数据结果：\n")
      print(res)      
      cat("------------------------\n")
      cat("--- 组间分析统计结果 ---\n")
      cat("分析方法：Kruskal-Wallis test", "\n")
      cat("组间自由度 (df)：", df_val, "\n")
      cat("卡方值 (Chi-squared)：", round(chi_sq, 3), "\n")
      cat("P 值：", format.pval(p_val, digits = 4, eps = 0.001), "\n")
      
      if(p_val < 0.05){
        cat(sprintf('结论：可以认为基因 [%s] 在不同分组间的表达水平分布显著不同（显著差异存在于中位数）。',
                    input$between_group_test_gene))
      } else {
        cat(sprintf('结论：目前没有足够证据表明基因 [%s] 在各组间的表达分布存在显著差异。',
                    input$between_group_test_gene))
      }
    }
  })
  
  
  ##事后检验server#####
  #加载事后检验基因
  observeEvent(result(),{
    req(result())
    updateSelectInput(inputId = 'post_hoc_gene',choices = unique(result()[[input$gene]]))
  })
  
  #渲染说明文本
  output$Explanation_text <- renderUI({
    tagList(
      
      
      # 参数检验部分
      tags$h4(tags$span(class = "label label-primary", "参数检验"), 
              tags$small(" (适用于正态且方差齐性数据/ANOVA显著后)")),
      tags$ul(class = "list-group",
              tags$li(class = "list-group-item", tags$strong("Tukey HSD:"), " 全组两两比较，学术最常用，严格控制假阳性。"),
              tags$li(class = "list-group-item", tags$strong("Dunnett:"), " 仅对比对照组 (Control vs Others)，统计效能最高。"),
              tags$li(class = "list-group-item", tags$strong("Bonferroni:"), " 最保守的校正，适用于比较次数较少的情况。"),
              tags$li(class = "list-group-item", tags$strong("LSD:"), " 最灵敏，但易产生假阳性，请谨慎引用。")
      ),
      
      # 非参数检验部分
      tags$h4(tags$span(class = "label label-warning", "非参数检验"), 
              tags$small(" (适用于偏态或方差不齐数据/Kruskal-Wallis显著后)")),
      tags$ul(class = "list-group",
              tags$li(class = "list-group-item", tags$strong("Dunn's Test:"), " 非参数两两比较的金标准，推荐首选。"),
              tags$li(class = "list-group-item", tags$strong("Wilcoxon Pairwise:"), " 基于秩次的配对比较，适合灵活的自定义校正。"),
              tags$li(class = "list-group-item", tags$strong("Nemenyi Test:"), " 适用于多组样本量均衡时的非参数比较。")
      ),
      
      # 注意事项提示框
      tags$div(class = "alert alert-info", style = "margin-top: 10px; padding: 10px;",
               tags$i(class = "fa fa-info-circle"),
               tags$strong(" 注意："), 
               "请务必先参考“方差齐性检验”结果。若方差不齐（P < 0.05），强烈建议优先选择非参数检验方法。"
      )
    )
  })
  
  #渲染每个选择方法的UI
  output$dynamic_post_hoc_results <- renderUI({
    # 只有点击了“开始检验”按钮且有选择方法时才执行
    req(input$start_post_hoc)
    methods <- input$post_hoc_methods
    req(length(methods) > 0)
    
    # 为每个选中的方法创建一个 UI 块
    tagList(
      lapply(methods, function(m) {
        # 为每种方法定义不同的颜色或标题
        box_color <- ifelse(m %in% c("Tukey", "Dunnett", "Bonferroni", "LSD"), "primary", "warning")
        
        box(
          title = paste("分析结果:", m),
          status = box_color,
          solidHeader = TRUE,
          width = 6, # 占据整行，也可以设为 6 并排显示
          collapsible = TRUE,
          height = "500px",
          # 这里的 ID 需要与后面的 renderPrint/renderDataTable 对应
          # 使用 prefix + method name 的方式命名
          tags$div(
            style = "height: 400px; overflow-y: auto;",
            dataTableOutput(outputId = paste0("res_out_", m))
          ),
          
        )
      })
    )
  })
  
  #事后检验
  post_hoc_results_list <- eventReactive(input$start_post_hoc,{
    req(result(),
        input$group,
        input$gene,
        input$reference_group,
        input$post_hoc_methods,
        input$post_hoc_gene)
    # 锁定当前计算的基因名
    target_gene <- input$post_hoc_gene
    
    #过滤数据
    sub_df <- result()[result()[[input$gene]] == target_gene, ]
    sub_df[[input$group]] <- as.factor(sub_df[[input$group]])
    
    # 存储所有选中方法的结果
    all_res <- list()
    
    for( m in input$post_hoc_methods){
      
      res <- tryCatch({
        # 构建统一的公式
        formula_obj <- as.formula(sprintf("`result` ~ `%s`", input$group))
        
        # 根据方法字符串调用函数
        result_obj <- switch (m,
          'Tukey'      = PMCMRplus::tukeyTest(formula_obj,data = sub_df),
          'LSD'        = PMCMRplus::lsdTest(formula_obj, data = sub_df),
          'Dunnett'    = PMCMRplus::dunnettTest(formula_obj, data = sub_df, ref = input$reference_group),
          "Dunn"       = PMCMRplus::kwAllPairsDunnTest(formula_obj,data = sub_df,p.adjust.method = "holm"),
          'Nemenyi'    = PMCMRplus::kwAllPairsNemenyiTest(formula_obj,data = sub_df),
          "Wilcoxon"   = pairwise.wilcox.test(sub_df[['result']], sub_df[[input$group]], p.adjust.method = "BH"),
          "Bonferroni" = pairwise.t.test(sub_df[['result']],sub_df[[input$group]],p.adjust.method = "bonferroni")
        )
        result_obj
      },error = function(e){ paste('统计计算失败：',e$message)})
      all_res[[m]] <- res
    }
    return(list(
      data_results = all_res,     # 具体的统计结果
      calc_gene = target_gene))    # 计算时对应的基因（属性）
  })
  
  # 动态分发渲染事后检验结果
  observe({
    # 确保计算已经完成
    req(post_hoc_results_list())
    full_obj <- post_hoc_results_list()
    
    # 提取计算时的基因属性
    current_gene_attr <- full_obj$calc_gene
    results_all <- full_obj$data_results
    
    methods_selected <- names(results_all)
    
    # 遍历选中的每一个方法
    for (m in methods_selected) {
      # 必须使用 local 锁定当前的 m，否则循环结束后所有输出都会指向最后一个方法
      local({
        method_name <- m
        # 注意：这里我们依然使用渲染时的基因属性，确保数据与标题一致
        gene_label <- current_gene_attr 
        
        output[[paste0("res_out_", method_name)]] <- renderDataTable({
         
          res_raw <- results_all[[method_name]]
          
          # 判断：如果计算失败（返回的是字符串），则显示错误提示
          if (is.character(res_raw)) {
            return(datatable(data.frame(错误 = res_raw)))
          }
          
          clean_df <- format_posthoc_res(res_raw)
          
          # 使用 DT 渲染美化表格
          DT::datatable(clean_df, 
                        caption = paste("基因:", full_obj$calc_gene, " 的 ", method_name, " 检验结果"),
                        options = list(dom = 't', 
                                       pageLength = -1,
                                       columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                       paging = FALSE,
                                       initComplete = JS(
                                         "function(settings, json) {",
                                         "$(this.api().table().header()).css({'background-color': '#f8f9fa', 'color': '#333', 'text-align': 'center'});",
                                         "}"
                                       )
                                       ), 
                        selection = 'none',
                        rownames = FALSE) %>%
            DT::formatStyle('显著性', 
                            color = DT::styleEqual(c("***", "**", "*"), c('red', 'red', 'orange')))
          
          
        })
      })
    }
  })
  
  #下载事后检验结果
  output$download_all_posthoc <- downloadHandler(
    filename = function(){
      gene_name <- post_hoc_results_list()$calc_gene
      paste0('事后检验_',gene_name,'_',format(Sys.time(),'%Y%m%d-%H%M'),'.xlsx')
    },
    content = function(file){
      req(post_hoc_results_list())
      full_obj <- post_hoc_results_list()
      results_all <- full_obj$data_results
      
      # 创建一个新的 Excel 工作簿
      wb <- createWorkbook()
      
      for( m in names(results_all)){
        res_raw <- results_all[[m]]
        
        if(!is.character(res_raw)){
          clean_df <- format_posthoc_res(res_raw)
          # 添加工作表
          addWorksheet(wb,sheetName = m)
          
          #写入数据
          writeData(wb,sheet = m,x=clean_df)
          
          #可选：美化 Excel 表头样式
          headerStyle <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1")
          addStyle(wb, sheet = m, headerStyle, rows = 1, cols = 1:ncol(clean_df), gridExpand = TRUE)
        }
      }
      saveWorkbook(wb,file,overwrite = T)
    }
  )
  
  
  
  
  
  
  ##绘图server#####
  
  #更新基因选项
  observeEvent(list(result(), input$gene),{
    req(result(),input$gene)
    updateSelectInput(inputId = 'plot_gene',choices = unique(result()[[input$gene]]))
  })
  
  
  # 根据加载的数据动态生成排序选择器
  output$group_order_ui <- renderUI({
    req(result(), input$group)
    # 获取当前数据中所有的分组名称
    groups <- unique(result()[[input$group]])
    
    # 使用  pickerInput 允许用户拖拽或按顺序点击
    selectInput(
      inputId = "group_order",
      label = "自定义分组排序 (按顺序点击选择):",
      choices = groups,
      selected = groups, # 默认顺序
      multiple = TRUE,
      width = "100%"
    )
  })
  
  
  ##生成ggplot2对象
  plot_obj <- eventReactive(input$render_plot,{
    req(result(),input$plot_gene)
    sub_df <- result()[result()[[input$gene]] == input$plot_gene, ]
    
    if (!is.null(input$group_order) && length(input$group_order) > 0) {
      # 只针对出现在用户选择列表里的组进行排序
      sub_df[[input$group]] <- factor(
        sub_df[[input$group]], 
        levels = input$group_order
      )}
    
    draw_plot(
      data        = sub_df,
      x_var       = input$group,
      y_var       = "result", 
      plot_type   = input$plot_type,
      palette     = input$plot_palette,
      bin_width   = input$bin_width,      
      xlabs       = input$xlab,
      ylabs       = input$ylab,
      font_family = input$font_family,
      font_size   = input$font_size,
      alpha       = if(input$fill_transparent) input$plot_alpha else 1,
      show_points = input$show_points,
      show_points_color = input$show_points_color
    )
  })
  ##渲染
  output$main_plot <- renderPlot({
    plot_obj()
  })
  
  #图片下载
  output$download_plot <- downloadHandler(
    filename = function(){
      paste0('qPCR_plot_',input$plot_gene,format(Sys.time(),"%Y%m%d_%H%M"),'.',input$file_type)
    },
    content = function(file){
      req(plot_obj())
      
      p <- plot_obj()
      
      
      ggsave(
        filename = file,
        plot     = p,
        device   = input$file_type,
        width    = input$plot_width,
        height   = input$plot_height,
        dpi      = input$plot_res,
        units    = 'in'
      )
    }
  )
}

shinyApp(ui, server)
























