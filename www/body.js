$(document).ready(function(){
  
  $('body').addClass("control-sidebar-open");
  
  var content_width = $(".content").css('width');
  var rightbar_width = $("#controlbar aside").css('width');
  var size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
  var size_plot = (size_table / 2) - 10;
  
  resize_box(size_table, size_plot);
  
  pp_parent_w = $("#pca_plot_box").parent('div').width();
  $("#pca_plot_box").parent('div').css("height",pp_parent_w*0.8);
  cm_parent_w = $("#correlation_matrix_box").parent('div').width();
  $("#correlation_matrix_box").parent('div').css("height",cm_parent_w*0.8);
  hm_parent_w = $("#heatmap_box").parent('div').width();
  $("#heatmap_box").parent('div').css("height",hm_parent_w*0.8);

  $( window ).resize( function() {
    var window_width = $(window).width();
    if(window_width >= 2000){
      content_width = $(".content").css('width');
      rightbar_width = $("#controlbar aside").css('width');
      
      size_table = window_width - (parseInt(rightbar_width)*1.9);
      size_plot = (size_table / 2)-10;
      resize_box(size_table, size_plot);
    }else{
      content_width = $(".content").css - 20;
      rightbar_width = $("#controlbar aside").css('width');
      content_width = window_width - (parseInt(rightbar_width)*1.9);
      
      size_table = parseInt(content_width);
      size_plot = (size_table / 2) - 10;
      resize_box(size_table, size_plot);
    }
  });
  

  $("#gobp_gsa_box .load-container").hide();
  $("#gocc_gsa_box .load-container").hide();
  $("#gomf_gsa_box .load-container").hide();
  $("#kegg_gsa_box .load-container").hide();
  
  $("#gsa_btn").click(function(){
    $("#gobp_gsa_box .load-container").show();
    $("#gocc_gsa_box .load-container").show();
    $("#gomf_gsa_box .load-container").show();
    $("#kegg_gsa_box .load-container").show(); 
  });

  $("#gobp_gsea_box .load-container").hide();
  $("#gocc_gsea_box .load-container").hide();
  $("#gomf_gsea_box .load-container").hide();
  $("#kegg_gsea_box .load-container").hide();
  
  $("#gsea_btn").click(function(){
    $("#gobp_gsea_box .load-container").show();
    $("#gocc_gsea_box .load-container").show();
    $("#gomf_gsea_box .load-container").show();
    $("#kegg_gsea_box .load-container").show();
  });

  
  /*$("#ppi_box .load-container").hide();
  $("#ppi_btn").click(function(){
     $("#ppi_box .load-container").show();
   });
  
  $("#dea_btn").click(function(){
    $(window).bind('load',function(){
      content_width = $(".content").css('width');
      rightbar_width = $("#controlbar aside").css('width');
      size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
      size_plot = (size_table / 2) - 10;
    
      resize_box(size_table, size_plot);
    });
  });*/
  
  var dimension = [0, 0];
  $("#zoom_pathway_btn").click(function(){
    dimension[0] = 1344;
    dimension[1] = 756;
    Shiny.onInputChange("dimension", dimension);
  });
 
  $(".nav.navbar-nav a").click(function(){
    content_width = $(".content").css('width');
    rightbar_width = $("#controlbar aside").css('width');
    size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
    size_plot = (size_table / 2) - 10;
    
    resize_box(size_table, size_plot);
  });
});

function resize_box(size_table, size_plot){
  bodyClass = $('body').attr('class');
    if(bodyClass == 'skin-blue sidebar-mini control-sidebar-open'){
      
      dt_parent1 = $("#data_table").parent('div');
      dt_parent2.css('width',size_table);
      dt_parent2 = dt_parent1.parent('div');
      
      pt_parent1 = $("#plot_tabBox").parent('div');
      pt_parent1.css('height',$("#plot_tabBox").css('height'));
      pt_parent2 = pt_parent1.parent('div');
      pt_parent2.css('width',size_table);
      
      //pca_plot_box
      pp_parent1 = $("#pca_plot_box").parent('div');
      //pp_parent1.css('height',$("#pca_plot_box").css('height'));
      pp_parent2 = pp_parent1.parent('div');
      pp_parent2.css('width',size_plot);
      
      //volcano_box
      vp_parent1 = $("#volcano_box").parent('div');
      var vp_height = $("#volcano_box").css('height')+$("#volcano_info").css('height')+$(".download_volcano").css('height');
      vp_parent1.css('height',vp_height);
      vp_parent2 = vp_parent1.parent('div');
      vp_parent2.css('width',size_plot);
      
      //correlation_matrix_box
      cm_parent1 = $("#correlation_matrix_box").parent('div');
      //cm_parent1.css('height',$("#correlation_matrix_box").css('height'));
      cm_parent2 = cm_parent1.parent('div');
      cm_parent2.css('width',size_plot);
      
      //heatmap_box
      hm_parent1 = $("#heatmap_box").parent('div');
      //hm_parent1.css('height',$("#heatmap_box").css('height'));
      hm_parent2 = hm_parent1.parent('div');
      hm_parent2.css('width',size_plot);
      
    } else{
       //data_table
      dt_parent1 = $("#data_table").parent('div');
      dt_parent2 = dt_parent1.parent('div');
      dt_parent2.css('width','');
      
      pt_parent1 = $("#plot_tabBox").parent('div');
      pt_parent2 = pt_parent1.parent('div');
      pt_parent2.css('width','');
      
      //pca_plot_box
      pp_parent1 = $("#pca_plot_box").parent('div');
      pp_parent2 = pp_parent1.parent('div');
      pp_parent2.css('width','');
      
      //volcano_box
      vp_parent1 = $("#volcano_box").parent('div');
      vp_parent2 = vp_parent1.parent('div');
      vp_parent2.css('width','');
      
      //correlation_matrix_box
      cm_parent1 = $("#correlation_matrix_box").parent('div');
      cm_parent2 = cm_parent1.parent('div');
      cm_parent2.css('width','');
      
      //heatmap_box
      hm_parent1 = $("#heatmap_box").parent('div');
      hm_parent2 = hm_parent1.parent('div');
      hm_parent2.css('width','');
    }
}