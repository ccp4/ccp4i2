require.config({
               //baseUrl: 'http://127.0.0.1:43434/report_files/0.1.0',
               paths: {
               // the left side is the module ID,
               // the right side is the path to
               // the jQuery file, relative to baseUrl.
               // Also, the path should NOT include
               // the '.js' file extension. This example
               // is using jQuery 1.9.0 located at
               // js/lib/jquery-1.9.0.js, relative to
               // the HTML page.
               jquery: 'dist/libs/jquery',
               mosflm: 'mosflm',
               flot: 'flot/jquery.flot',
               CCP4i2Widgets: 'CCP4i2Widgets',
               jspimple: 'jspimple',
               xreport:'xreport',
               ProjectFormatter:'ProjectFormatter',
               jstree:'jstree.min'
               }
               });

require(['jquery','xreport','jstree','ProjectFormatter','jspimple','jstables'],function($,xreport,jstree,jspimple,jstables){
        $( document ).ready(function(){onLoad();})
        })