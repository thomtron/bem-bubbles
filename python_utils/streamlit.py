import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mpld3
import streamlit.components.v1 as components
from mpld3 import plugins




css = """
table
{
  border-collapse: collapse;
}
th
{
  color: #ffffff;
  background-color: #000000;
}
td
{
  background-color: #cccccc;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 1px solid black;
  text-align: right;
}
"""

#create your figure and get the figure object returned
fig = plt.figure() 
ax = fig.add_subplot()
x = np.linspace(0,np.pi,100)
ax.plot(x,np.sin(x),'-k.',markersize=20) 


for line in ax.get_lines():
    xy_data = line.get_xydata()
    labels = []
    for x, y in xy_data:
        # Create a label for each point with the x and y coords
        html_label = f'<table border="1" class="dataframe"> <thead> <tr style="text-align: right;"> </thead> <tbody> <tr> <th>x</th> <td>{x}</td> </tr> <tr> <th>y</th> <td>{y}</td> </tr> </tbody> </table>'
        labels.append(html_label)

    # Create the tooltip with the labels (x and y coords) and attach it to each line with the css specified
    tooltip = plugins.PointHTMLTooltip(line, labels, css=css)
    # Since this is a separate plugin, you have to connect it
    plugins.connect(fig, tooltip)
    #plugins.connect(fig,ClickInfo(line))


fig_html = mpld3.fig_to_html(fig)
components.html(fig_html, height=600)
