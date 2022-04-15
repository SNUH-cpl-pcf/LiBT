$(document).ready(function(){
  
  var content_width = $(".content").css('width');
  var rightbar_width = $("#controlbarId").css('width');
  var size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
  var size_plot = (size_table / 2) - 10;
  
  resize_box(size_table, size_plot);

  $( window ).resize( function() {
    var window_width = $(window).width();
    if(window_width >= 2000){
      content_width = $(".content").css('width');
      rightbar_width = $("#controlbarId").css('width');
      
      size_table = window_width - (parseInt(rightbar_width)*1.9);
      size_plot = (size_table / 2)-10;
      resize_box(size_table, size_plot);
    }else{
      content_width = $(".content").css - 20;
      rightbar_width = $("#controlbarId").css('width');
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
    rightbar_width = $("#controlbarId").css('width');
    size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
    size_plot = (size_table / 2) - 10;
    
    resize_box(size_table, size_plot);
  });
  
  $(".sidebar-menu li a").click(function(){
    content_width = $(".content").css('width');
    rightbar_width = $("#controlbarId").css('width');
    size_table = parseInt(content_width) - parseInt(rightbar_width) - 20;
    size_plot = (size_table / 2) - 10;
    
    resize_box(size_table, size_plot);
  });
});

function resize_box(size_table, size_plot){
  var asideClass = $('#controlbarId').attr('class');
    if(asideClass == 'control-sidebar control-sidebar-dark shiny-bound-input control-sidebar-open' || asideClass == 'control-sidebar control-sidebar-dark control-sidebar-open shiny-bound-input'){
      
      var dt_parent1 = $("#data_table").parent('div');
      var dt_parent2 = dt_parent1.parent('div');
      dt_parent2.css('width',size_table);
      
      var pt_parent1 = $("#plot_tabBox").parent('div');
      pt_parent1.css('height',$("#plot_tabBox").css('height'));
      var pt_parent2 = pt_parent1.parent('div');
      pt_parent2.css('width',size_table);
      
    } else{
       //data_table
      dt_parent1 = $("#data_table").parent('div');
      dt_parent2 = dt_parent1.parent('div');
      dt_parent2.css('width','');
      
      pt_parent1 = $("#plot_tabBox").parent('div');
      pt_parent2 = pt_parent1.parent('div');
      pt_parent2.css('width','');
     
    }
}
