import os, sys
import numpy as np
import params
import pandas as pd

from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource
from bokeh.models.widgets import CheckboxButtonGroup, Slider
from bokeh.layouts import gridplot, column, row

samples = sys.argv[1]
output_file('final_interactive_'+samples+'.html')
start = 2.7 if samples == 'inner' else 1.45
end = 3 if samples == 'inner' else 1.62
y_range = (-1.5, 1.5)
x_range = (-1, 6.5)
colors = ['blue', 'darkorange', 'forestgreen', 'crimson']
colors_diamonds = ['royalblue', 'khaki', 'yellowgreen', 'red']

def create_final_box_plot(method, region):
    df, p1 = ([] for _ in range(2))
    masks = [3,4,5,6]
    for imask in masks:
        with np.load('numpy_files/arrays_reg'+region+'_'+str(imask)+samples+'_'+method+'.npz') as f:
            df.append(pd.DataFrame(f['response_eta']).transpose())
            df[-1].columns = ['response', 'eta']

    #data sources
    data_limit = 64800
    source1 = ColumnDataSource(data=dict(x_mask3=df[0].eta[:data_limit], x_mask4=df[1].eta[:data_limit], x_mask5=df[2].eta[:data_limit], x_mask6=df[3].eta[:data_limit], 
                                         y_mask3=df[0].response[:data_limit], y_mask4=df[1].response[:data_limit], y_mask5=df[2].response[:data_limit], y_mask6=df[3].response[:data_limit],
                                         c_mask3=[colors[0]]*data_limit, c_mask4=[colors[1]]*data_limit, c_mask5=[colors[2]]*data_limit, c_mask6=[colors[3]]*data_limit))
    source1_v2 = ColumnDataSource(data=dict(x_mask3=df[0].eta[:data_limit], x_mask4=df[1].eta[:data_limit], x_mask5=df[2].eta[:data_limit], x_mask6=df[3].eta[:data_limit],
                                            y_mask3=df[0].response[:data_limit], y_mask4=df[1].response[:data_limit], y_mask5=df[2].response[:data_limit], y_mask6=df[3].response[:data_limit],
                                            c_mask3=[colors[0]]*data_limit, c_mask4=[colors[1]]*data_limit, c_mask5=[colors[2]]*data_limit, c_mask6=[colors[3]]*data_limit))
    source_vlines = ColumnDataSource(data=dict(etamin=[start], etamax=[end]))
    source_outliers = ColumnDataSource(data=dict(x_mask3=[1]*data_limit, x_mask4=[2]*data_limit, x_mask5=[3]*data_limit, x_mask6=[4]*data_limit, 
                                                 y_mask3=df[0].response[:data_limit]+100., y_mask4=df[1].response[:data_limit]+100., y_mask5=df[2].response[:data_limit]+100., y_mask6=df[3].response[:data_limit]+100.))
    source2 = ColumnDataSource(data=dict(x_mask3=[1],x_mask4=[2],x_mask5=[3],x_mask6=[4], 
                                         y1_mask3=[0.],y1_mask4=[0.],y1_mask5=[0.],y1_mask6=[0.],
                                         y2_mask3=[0.], y2_mask4=[0.], y2_mask5=[0.], y2_mask6=[0.],
                                         height1_mask3=[.6],height1_mask4=[.6],height1_mask5=[.6],height1_mask6=[.6],
                                         height2_mask3=[.6],height2_mask4=[.6],height2_mask5=[.6],height2_mask6=[.6],
                                         q1_mask3=[-.25],q1_mask4=[-.25],q1_mask5=[-.25],q1_mask6=[-.25], 
                                         q3_mask3=[-.25],q3_mask4=[-.25],q3_mask5=[-.25],q3_mask6=[-.25], 
                                         upper_mask3=[-.6],upper_mask4=[-.6],upper_mask5=[-.6],upper_mask6=[-.6],
                                         lower_mask3=[.6], lower_mask4=[.6], lower_mask5=[.6], lower_mask6=[.6]))

    #1st figure
    plot_options = dict(plot_height=500, plot_width=350, y_range=y_range, tools="wheel_zoom,box_zoom,box_select,pan,reset", output_backend="webgl")
    for imask in masks:
        p1.append(figure(title="Eta vs. Response          Mask "+str(imask), **plot_options))
        p1[-1].circle('x_mask'+str(imask), 'y_mask'+str(imask), color='c_mask'+str(imask), size=2, alpha=0.4, source=source1)
        p1[-1].segment(x0='etamin', y0=y_range[0], x1='etamin', y1=y_range[1], source=source_vlines, line_color='black', line_width=2)
        p1[-1].segment(x0='etamax', y0=y_range[0], x1='etamax', y1=y_range[1], source=source_vlines, line_color='black', line_width=2)

    #2nd figure
    plot2_options = dict(plot_height=300, plot_width=500, x_range=x_range, y_range=(y_range[0]-0.3,y_range[1]+0.3), tools="wheel_zoom,box_zoom,pan,reset", output_backend="webgl")
    radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
    methoddict = dict({'nocorr':'No correction', 'corr_ed':'Shower leakage', 'corr_fineeta':'Brute force'})
    p2 = figure(title="Signal integration radius: "+radius[region]+'        Method: '+methoddict[method], **plot2_options)
    p2.xaxis.visible = False
    box_options = dict(line_color='black', line_width=2, source=source2)
    for imask in masks:
        #boxes
        p2.rect(x='x_mask'+str(imask), y='y1_mask'+str(imask), width=0.9, height='height1_mask'+str(imask), color=colors[imask-masks[0]], legend='Mask '+str(imask), **box_options)
        p2.rect(x='x_mask'+str(imask), y='y2_mask'+str(imask), width=0.9, height='height2_mask'+str(imask), color=colors[imask-masks[0]], **box_options)
        #segments
        p2.segment(x0='x_mask'+str(imask), y0='q3_mask'+str(imask), x1='x_mask'+str(imask), y1='upper_mask'+str(imask), **box_options)
        p2.segment(x0='x_mask'+str(imask), y0='q1_mask'+str(imask), x1='x_mask'+str(imask), y1='lower_mask'+str(imask), **box_options)
        #whiskers
        p2.rect(x='x_mask'+str(imask), y='upper_mask'+str(imask), width=0.1, height=0.005, **box_options)
        p2.rect(x='x_mask'+str(imask), y='lower_mask'+str(imask), width=0.1, height=0.005, **box_options)
        #outliers
        p2.diamond(x='x_mask'+str(imask), y='y_mask'+str(imask), line_color=colors_diamonds[imask-masks[0]], size=4, source=source_outliers)
        #dashed line
        p2.segment(x0=x_range[0]+0.5, y0=0., x1=x_range[1]-0.5, y1=0., line_color='black', line_dash='dashed')

    #Dynamic behaviour (the model that triggers the callback is called 'cb_obj')
    callback = CustomJS(args=dict(source1=source1, source1_v2=source1_v2, source2=source2, 
                                  source_vlines=source_vlines, source_outliers=source_outliers), code="""
            var etamin = cb_obj.value;
            var delta_eta = 0.03;
            var data1 = source1.data;
            var data1_v2 = source1_v2.data;
            var data2 = source2.data;
            var data_lines = source_vlines.data;
            var data_out = source_outliers.data;
            data_lines['etamin'][0] = etamin;
            data_lines['etamax'][0] = etamin+delta_eta;

            data1_v2['x_mask3'] = [];
            data1_v2['x_mask4'] = [];
            data1_v2['x_mask5'] = [];
            data1_v2['x_mask6'] = [];
            data1_v2['y_mask3'] = [];
            data1_v2['y_mask4'] = [];
            data1_v2['y_mask5'] = [];
            data1_v2['y_mask6'] = [];
            for (var i = 0; i < data1['x_mask3'].length; i++) {
                if (data1['x_mask3'][i] > etamin && data1['x_mask3'][i] <= etamin+delta_eta) {
                   data1_v2['x_mask3'].push(data1['x_mask3'][i]);
                   data1_v2['y_mask3'].push(data1['y_mask3'][i]);
                }
            }
            for (var i = 0; i < data1['x_mask4'].length; i++) {
                if (data1['x_mask4'][i] > etamin && data1['x_mask4'][i] <= etamin+delta_eta) {
                   data1_v2['x_mask4'].push(data1['x_mask4'][i]);
                   data1_v2['y_mask4'].push(data1['y_mask4'][i]);
                }
            }
            for (var i = 0; i < data1['x_mask5'].length; i++) {
                if (data1['x_mask5'][i] > etamin && data1['x_mask5'][i] <= etamin+delta_eta) {
                   data1_v2['x_mask5'].push(data1['x_mask5'][i]);
                   data1_v2['y_mask5'].push(data1['y_mask5'][i]);
                }
            }
            for (var i = 0; i < data1['x_mask6'].length; i++) {
                if (data1['x_mask6'][i] > etamin && data1['x_mask6'][i] <= etamin+delta_eta) {
                   data1_v2['x_mask6'].push(data1['x_mask6'][i]);
                   data1_v2['y_mask6'].push(data1['y_mask6'][i]);
                }
            }

           //define quantiles and related quantities
           var data_tmp_mask3 = data1_v2['y_mask3'];
           data_tmp_mask3.sort( function(a,b) {return a - b;} );
           var l1_3 = Math.floor( (data_tmp_mask3.length) * 0.25);
           var l2_3 = Math.floor( (data_tmp_mask3.length) * 0.50);
           var l3_3 = Math.floor( (data_tmp_mask3.length) * 0.75);

           var data_tmp_mask4 = data1_v2['y_mask4'];
           data_tmp_mask4.sort( function(a,b) {return a - b;} );
           var l1_4 = Math.floor( (data_tmp_mask4.length) * 0.25);
           var l2_4 = Math.floor( (data_tmp_mask4.length) * 0.50);
           var l3_4 = Math.floor( (data_tmp_mask4.length) * 0.75);

           var data_tmp_mask5 = data1_v2['y_mask5'];
           data_tmp_mask5.sort( function(a,b) {return a - b;} );
           var l1_5 = Math.floor( (data_tmp_mask5.length) * 0.25);
           var l2_5 = Math.floor( (data_tmp_mask5.length) * 0.50);
           var l3_5 = Math.floor( (data_tmp_mask5.length) * 0.75);

           var data_tmp_mask6 = data1_v2['y_mask6'];
           data_tmp_mask6.sort( function(a,b) {return a - b;} );
           var l1_6 = Math.floor( (data_tmp_mask6.length) * 0.25);
           var l2_6 = Math.floor( (data_tmp_mask6.length) * 0.50);
           var l3_6 = Math.floor( (data_tmp_mask6.length) * 0.75);

           var upper_3;
           var lower_3;
           var upper_4;
           var lower_4;
           var upper_5;
           var lower_5;
           var upper_6;
           var lower_6;

           if(data_tmp_mask3.length % 2) 
           {
              var q1 = data_tmp_mask3[l1_3];
              var q2 = data_tmp_mask3[l2_3];
              var q3 = data_tmp_mask3[l3_3];
              data2['y1_mask3'][0] = (q2+q1)/2;
              data2['y2_mask3'][0] = (q3+q2)/2;
              data2['height1_mask3'][0] = q2-q1;
              data2['height2_mask3'][0] = q3-q2;
              data2['q1_mask3'][0] = q1;
              data2['q3_mask3'][0] = q3;
              upper_3 = q3+1.5*(q3-q1);
              lower_3 = q1-1.5*(q3-q1);
              data2['upper_mask3'][0] = upper_3;
              data2['lower_mask3'][0] = lower_3;
           }
           else
           {
              var q1 = (data_tmp_mask3[l1_3-1]+data_tmp_mask3[l1_3])/2.0;
              var q2 = (data_tmp_mask3[l2_3-1]+data_tmp_mask3[l2_3])/2.0;
              var q3 = (data_tmp_mask3[l3_3-1]+data_tmp_mask3[l3_3])/2.0;
              data2['y1_mask3'][0] = (q2+q1)/2;
              data2['y2_mask3'][0] = (q3+q2)/2;
              data2['height1_mask3'][0] = q2-q1;
              data2['height2_mask3'][0] = q3-q2;
              data2['q1_mask3'][0] = q1;
              data2['q3_mask3'][0] = q3;
              upper_3 = q3+1.5*(q3-q1);
              lower_3 = q1-1.5*(q3-q1);
              data2['upper_mask3'][0] = upper_3;
              data2['lower_mask3'][0] = lower_3;
           }

           if(data_tmp_mask4.length % 2) 
           {
              var q1 = data_tmp_mask4[l1_4];
              var q2 = data_tmp_mask4[l2_4];
              var q3 = data_tmp_mask4[l3_4];
              data2['y1_mask4'][0] = (q2+q1)/2;
              data2['y2_mask4'][0] = (q3+q2)/2;
              data2['height1_mask4'][0] = q2-q1;
              data2['height2_mask4'][0] = q3-q2;
              data2['q1_mask4'][0] = q1;
              data2['q3_mask4'][0] = q3;
              upper_4 = q3+1.5*(q3-q1);
              lower_4 = q1-1.5*(q3-q1);
              data2['upper_mask4'][0] = upper_4;
              data2['lower_mask4'][0] = lower_4;
           }
           else
           {
              var q1 = (data_tmp_mask4[l1_4-1]+data_tmp_mask4[l1_4])/2.0;
              var q2 = (data_tmp_mask4[l2_4-1]+data_tmp_mask4[l2_4])/2.0;
              var q3 = (data_tmp_mask4[l3_4-1]+data_tmp_mask4[l3_4])/2.0;
              data2['y1_mask4'][0] = (q2+q1)/2;
              data2['y2_mask4'][0] = (q3+q2)/2;
              data2['height1_mask4'][0] = q2-q1;
              data2['height2_mask4'][0] = q3-q2;
              data2['q1_mask4'][0] = q1;
              data2['q3_mask4'][0] = q3;
              upper_4 = q3+1.5*(q3-q1);
              lower_4 = q1-1.5*(q3-q1);
              data2['upper_mask4'][0] = upper_4;
              data2['lower_mask4'][0] = lower_4;
           }

           if(data_tmp_mask5.length % 2) 
           {
              var q1 = data_tmp_mask5[l1_5];
              var q2 = data_tmp_mask5[l2_5];
              var q3 = data_tmp_mask5[l3_5];
              data2['y1_mask5'][0] = (q2+q1)/2;
              data2['y2_mask5'][0] = (q3+q2)/2;
              data2['height1_mask5'][0] = q2-q1;
              data2['height2_mask5'][0] = q3-q2;
              data2['q1_mask5'][0] = q1;
              data2['q3_mask5'][0] = q3;
              upper_5 = q3+1.5*(q3-q1);
              lower_5 = q1-1.5*(q3-q1);
              data2['upper_mask5'][0] = upper_5;
              data2['lower_mask5'][0] = lower_5;
            }
            else
            {
              q1 = (data_tmp_mask5[l1_5-1]+data_tmp_mask5[l1_5])/2.0;
              q2 = (data_tmp_mask5[l2_5-1]+data_tmp_mask5[l2_5])/2.0;
              q3 = (data_tmp_mask5[l3_5-1]+data_tmp_mask5[l3_5])/2.0;
              data2['y1_mask5'][0] = (q2+q1)/2;
              data2['y2_mask5'][0] = (q3+q2)/2;
              data2['height1_mask5'][0] = q2-q1;
              data2['height2_mask5'][0] = q3-q2;
              data2['q1_mask5'][0] = q1;
              data2['q3_mask5'][0] = q3;
              upper_5 = q3+1.5*(q3-q1);
              lower_5 = q1-1.5*(q3-q1);
              data2['upper_mask5'][0] = upper_5;
              data2['lower_mask5'][0] = lower_5;
            }

           if(data_tmp_mask6.length % 2) 
           {
              var q1 = data_tmp_mask6[l1_6];
              var q2 = data_tmp_mask6[l2_6];
              var q3 = data_tmp_mask6[l3_6];
              data2['y1_mask6'][0] = (q2+q1)/2;
              data2['y2_mask6'][0] = (q3+q2)/2;
              data2['height1_mask6'][0] = q2-q1;
              data2['height2_mask6'][0] = q3-q2;
              data2['q1_mask6'][0] = q1;
              data2['q3_mask6'][0] = q3;
              upper_6 = q3+1.5*(q3-q1);
              lower_6 = q1-1.5*(q3-q1);
              data2['upper_mask6'][0] = upper_6;
              data2['lower_mask6'][0] = lower_6;
           }
           else
           {
              var q1 = (data_tmp_mask6[l1_6-1]+data_tmp_mask6[l1_6])/2.0;
              var q2 = (data_tmp_mask6[l2_6-1]+data_tmp_mask6[l2_6])/2.0;
              var q3 = (data_tmp_mask6[l3_6-1]+data_tmp_mask6[l3_6])/2.0;
              data2['y1_mask6'][0] = (q2+q1)/2;
              data2['y2_mask6'][0] = (q3+q2)/2;
              data2['height1_mask6'][0] = q2-q1;
              data2['height2_mask6'][0] = q3-q2;
              data2['q1_mask6'][0] = q1;
              data2['q3_mask6'][0] = q3;
              upper_6 = q3+1.5*(q3-q1);
              lower_6 = q1-1.5*(q3-q1);
              data2['upper_mask6'][0] = upper_6;
              data2['lower_mask6'][0] = lower_6;
           }

           //avoid putting the whiskers beyond the data points
           var currMax_3 = Math.max.apply(Math, data1_v2['y_mask3']);
           var currMin_3 = Math.min.apply(Math, data1_v2['y_mask3']);
           var currMax_4 = Math.max.apply(Math, data1_v2['y_mask4']);
           var currMin_4 = Math.min.apply(Math, data1_v2['y_mask4']);
           var currMax_5 = Math.max.apply(Math, data1_v2['y_mask5']);
           var currMin_5 = Math.min.apply(Math, data1_v2['y_mask5']);
           var currMax_6 = Math.max.apply(Math, data1_v2['y_mask6']);
           var currMin_6 = Math.min.apply(Math, data1_v2['y_mask6']);

           if(currMax_3 > upper_3)
              data2['upper_mask3'][0] = upper_3; 
           else
              data2['upper_mask3'][0] = currMax_3; 
           if(currMin_3 < lower_3)
              data2['lower_mask3'][0] = lower_3; 
           else
              data2['lower_mask3'][0] = currMin_3; 

           if(currMax_4 > upper_4)
              data2['upper_mask4'][0] = upper_4; 
           else
              data2['upper_mask4'][0] = currMax_4; 
           if(currMin_4 < lower_4)
              data2['lower_mask4'][0] = lower_4; 
           else
              data2['lower_mask4'][0] = currMin_4; 

           if(currMax_5 > upper_5)
              data2['upper_mask5'][0] = upper_5; 
           else
              data2['upper_mask5'][0] = currMax_5; 
            if(currMin_5 < lower_5)
              data2['lower_mask5'][0] = lower_5; 
           else
              data2['lower_mask5'][0] = currMin_5; 

           if(currMax_6 > upper_6)
               data2['upper_mask6'][0] = upper_6; 
           else
               data2['upper_mask6'][0] = currMax_6; 
           if(currMin_6 < lower_6)
               data2['lower_mask6'][0] = lower_6; 
           else
               data2['lower_mask6'][0] = currMin_6; 

           //place outliers
           data_out['y_mask3'] = [];
           data_out['y_mask4'] = [];
           data_out['y_mask5'] = [];
           data_out['y_mask6'] = [];
           for (var i = 0; i < data1_v2['x_mask3'].length; i++) {
             if (data1_v2['y_mask3'][i] > upper_3 || data1_v2['y_mask3'][i] < lower_3)
                  data_out['y_mask3'].push(data1_v2['y_mask3'][i]);
           }
           for (var i = 0; i < data1_v2['x_mask4'].length; i++) {
             if (data1_v2['y_mask4'][i] > upper_4 || data1_v2['y_mask4'][i] < lower_4)
                  data_out['y_mask4'].push(data1_v2['y_mask4'][i]);
           }
           for (var i = 0; i < data1_v2['x_mask5'].length; i++) {
             if (data1_v2['y_mask5'][i] > upper_5 || data1_v2['y_mask5'][i] < lower_5)
                  data_out['y_mask5'].push(data1_v2['y_mask5'][i]);
           }
           for (var i = 0; i < data1_v2['x_mask6'].length; i++) {
             if (data1_v2['y_mask6'][i] > upper_6 || data1_v2['y_mask6'][i] < lower_6)
                  data_out['y_mask6'].push(data1_v2['y_mask6'][i]);
           }

            //update sources
            source1.change.emit();
            source1_v2.change.emit();
            source2.change.emit();
            source_vlines.change.emit();
            source_outliers.change.emit();
        """)
    slider = Slider(start=start, end=end, value=start, step=.0001, title="Pseudorapidity (left bin edge)")
    slider.js_on_change('value', callback)
    return slider, p2


slider1a, box1a = create_final_box_plot('nocorr', '1')
slider2a, box2a = create_final_box_plot('corr_ed', '1')
slider3a, box3a = create_final_box_plot('corr_fineeta', '1')
slider1b, box1b = create_final_box_plot('nocorr', '2')
slider2b, box2b = create_final_box_plot('corr_ed', '2')
slider3b, box3b = create_final_box_plot('corr_fineeta', '2')
slider1c, box1c = create_final_box_plot('nocorr', '3')
slider2c, box2c = create_final_box_plot('corr_ed', '3')
slider3c, box3c = create_final_box_plot('corr_fineeta', '3')

sliders1 = [slider1a, slider2a, slider3a]
sliders2 = [slider1b, slider2b, slider3b]
sliders3 = [slider1c, slider2c, slider3c]
boxes1 = [box1a, box2a, box3a]
boxes2 = [box1b, box2b, box3b]
boxes3 = [box1c, box2c, box3c]

col1 = column(slider1a)
col2 = column(slider1b)
col3 = column(slider1c)
checkbox = CheckboxButtonGroup(labels=["No correction", "Shower leakage", "Brute force"], active=[])
checkbox_sr = CheckboxButtonGroup(labels=["1.3cm", "2.6cm", "5.3cm"], active=[])
checkbox_callback = CustomJS(args=dict(sliders1=sliders1, sliders2=sliders2, sliders3=sliders3, boxes1=boxes1, boxes2=boxes2, boxes3=boxes3, 
                                       col1=col1, col2=col2, col3=col3, checkbox=checkbox, checkbox_sr=checkbox_sr), code="""
                      const children_sr1 = [];
                      const children_sr2 = [];
                      const children_sr3 = [];
                      for (const i of checkbox.active) {
                         if (checkbox_sr.active.includes(0)) {
                            children_sr1.push(sliders1[i]);
                            children_sr1.push(boxes1[i]);
                         }
                         if (checkbox_sr.active.includes(1)) {
                            children_sr2.push(sliders2[i]);
                            children_sr2.push(boxes2[i]);
                         }
                         if (checkbox_sr.active.includes(2)) {
                            children_sr3.push(sliders3[i]);
                            children_sr3.push(boxes3[i]);
                         }
                      } 
                      col1.children = children_sr1;
                      col2.children = children_sr2;
                      col3.children = children_sr3;
                      """)

checkbox.js_on_change('active', checkbox_callback)
checkbox_sr.js_on_change('active', checkbox_callback)

layout = gridplot([[checkbox],[checkbox_sr], [col1, col2, col3]])
show(layout)
