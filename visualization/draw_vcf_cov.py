#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging
import gzip

# Third party imports
import argparse
import datetime
import pandas as pd

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: search for annotated vcf files and coverage files and outputs an interactive graph for every pair of files
INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 06 May 2020
REVISION: 

TODO:

================================================================
END_OF_HEADER
================================================================
"""

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def import_VCF41_to_pandas(vcf_file, sep='\t'):
    header_lines = 0
    if vcf_file.endswith(".gz"):
        compress = 'gzip'
        with gzip.open(vcf_file, 'rb') as f:
            first_line = f.readline().decode().strip()
            next_line = f.readline().decode().strip()
            while next_line.startswith("##"):
                header_lines = header_lines + 1
                next_line = f.readline().decode().strip()
    else:
        compress = None
        with open(vcf_file, 'r') as f:
            first_line = f.readline().strip()
            next_line = f.readline().strip()
            while next_line.startswith("##"):
                header_lines = header_lines + 1
                next_line = f.readline().strip()
    
    if first_line.endswith('VCFv4.1'):
        dataframe = pd.read_csv(vcf_file, compression=compress, sep=sep, skiprows=[header_lines], header=header_lines)

        sample = dataframe.columns[-1]
        dataframe.rename(columns={sample:'sample'}, inplace=True)
        
        for index, data_row in dataframe.iterrows():
            info_fields = re.findall(r';*([a-zA-Z]{1,20})=', data_row.INFO)
            info_values = re.findall(r'-?\d+\.?\d*e?[+-]?\d{0,2}', data_row.INFO)
            
            format_fields = data_row['FORMAT'].split(":")
            format_values = data_row['sample'].split(":")
                                    
            for ifield, ivalue in zip(info_fields,info_values):
                dataframe.loc[index,ifield] = ivalue
                
            for ffield, fvalue in zip(format_fields,format_values):
                dataframe.loc[index,ffield] = fvalue
                
            if 'ANN' in info_fields:
                annotation = dataframe.loc[index, 'INFO'].split('|')[1]
                dataframe.loc[index,'ANN'] = annotation
            else:
                dataframe.loc[index,'ANN'] = 'NONE'
                
        dataframe['FREQ'] = dataframe['FREQ'].str.replace('%', '')
                
        to_float = ['ADP', 'WT', 'HET', 'HOM', 'NC', 'GQ', 'SDP', 'DP',
                    'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']
            
        for column in dataframe.columns:
            if column in to_float:
                dataframe[column] = dataframe[column].astype(float)
                
        dataframe = dataframe[['REF','ALT','POS','FREQ', 'ANN']]
        
        
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe

def import_cov_to_pandas(cov_file):
    dataframe = pd.read_csv(cov_file, sep='\t', names=['REF', 'POS', 'DP'])
    dataframe = dataframe[['POS', 'DP']]
    return dataframe

def find_files(folder):
    vcf_files = []
    cov_files = []
    for root, _, files in os.walk(folder):
        for name in files:
            filename = os.path.join(root, name)
            if name.endswith('lowfreq.snpEff.vcf.gz'):
                vcf_files.append(filename)
            elif name.endswith('.cov'):
                cov_files.append(filename)
    vcf_files.sort()
    cov_files.sort()
    return vcf_files, cov_files

d3_template = """
    <!DOCTYPE html>
    <html>
        <head>
            <title>Allele_frequency_in_COVID19_samples</title>
            <style>
            html body {
                margin: 0 auto;
            }
            
            #center {
                display:flex;
                flex-direction:column;
                align-items:center;
            }

            text {
                font-family: sans-serif;
            }

            .tick text {
                font-size: 1em;
                fill: #635F5D;
            }

            .axisGrey .tick line{
                stroke: #C0C0BB;
            }
            
            .axisNone .tick line{
                stroke-opacity: 0;
            }

            .xAxisLabel{
                font-size: 2em;
                fill: #8E8883;
            }

            .yAxisLabel{
                font-size: 2em;
                fill: #8E8883;
            }

            .title{
                text-anchor: middle;
                font-size: 2em;
                fill: #635F5D;
            }
            
            div.tooltip {
                position: absolute;
                text-align: center;
                width: auto;
                height: auto;
                padding: 2px;
                font: 12px sans-serif;
                background: lightsteelblue;
                border: 0px;
                border-radius: 8px;
                pointer-events: none;
            }

            .line {
            fill: none;
            stroke: #6F257F;
            stroke-width: 2px;
            }
            
            .amplicon:hover {
                fill: brown;
            }
            
            .legend {
                font-size: 9px;
            }
            
            #apple{
                fill: black;
            }
        
            </style>
            <script src="https://d3js.org/d3.v5.min.js"></script>
        </head>
        <body>
            <div id="center">
                SVGGROUP
            </div>
                
            <script>

            const render = (selector, data, covdata, sample = "", svgwidth=1500, svgheight=500) => {
                
                beddata = [{'start': 30, 'end': 410, 'name': 'nCoV-2019_1'},
                        {'start': 320, 'end': 726, 'name': 'nCoV-2019_2'},
                        {'start': 642, 'end': 1028, 'name': 'nCoV-2019_3'},
                        {'start': 943, 'end': 1337, 'name': 'nCoV-2019_4'},
                        {'start': 1242, 'end': 1651, 'name': 'nCoV-2019_5'},
                        {'start': 1573, 'end': 1964, 'name': 'nCoV-2019_6'},
                        {'start': 1868, 'end': 2269, 'name': 'nCoV-2019_7_alt0'},
                        {'start': 2181, 'end': 2592, 'name': 'nCoV-2019_8'},
                        {'start': 2504, 'end': 2904, 'name': 'nCoV-2019_9_alt4'},
                        {'start': 2826, 'end': 3210, 'name': 'nCoV-2019_10'},
                        {'start': 3144, 'end': 3531, 'name': 'nCoV-2019_11'},
                        {'start': 3460, 'end': 3853, 'name': 'nCoV-2019_12'},
                        {'start': 3771, 'end': 4164, 'name': 'nCoV-2019_13'},
                        {'start': 4044, 'end': 4450, 'name': 'nCoV-2019_14_alt4'},
                        {'start': 4294, 'end': 4696, 'name': 'nCoV-2019_15'},
                        {'start': 4636, 'end': 5017, 'name': 'nCoV-2019_16'},
                        {'start': 4939, 'end': 5321, 'name': 'nCoV-2019_17'},
                        {'start': 5230, 'end': 5643, 'name': 'nCoV-2019_18'},
                        {'start': 5563, 'end': 5957, 'name': 'nCoV-2019_19'},
                        {'start': 5867, 'end': 6272, 'name': 'nCoV-2019_20'},
                        {'start': 6167, 'end': 6548, 'name': 'nCoV-2019_21'},
                        {'start': 6466, 'end': 6873, 'name': 'nCoV-2019_22'},
                        {'start': 6718, 'end': 7117, 'name': 'nCoV-2019_23'},
                        {'start': 7035, 'end': 7415, 'name': 'nCoV-2019_24'},
                        {'start': 7305, 'end': 7694, 'name': 'nCoV-2019_25'},
                        {'start': 7626, 'end': 8019, 'name': 'nCoV-2019_26'},
                        {'start': 7943, 'end': 8341, 'name': 'nCoV-2019_27'},
                        {'start': 8249, 'end': 8661, 'name': 'nCoV-2019_28'},
                        {'start': 8595, 'end': 8983, 'name': 'nCoV-2019_29'},
                        {'start': 8888, 'end': 9271, 'name': 'nCoV-2019_30'},
                        {'start': 9204, 'end': 9585, 'name': 'nCoV-2019_31'},
                        {'start': 9477, 'end': 9858, 'name': 'nCoV-2019_32'},
                        {'start': 9784, 'end': 10171, 'name': 'nCoV-2019_33'},
                        {'start': 10076, 'end': 10459, 'name': 'nCoV-2019_34'},
                        {'start': 10362, 'end': 10763, 'name': 'nCoV-2019_35'},
                        {'start': 10666, 'end': 11074, 'name': 'nCoV-2019_36'},
                        {'start': 10999, 'end': 11394, 'name': 'nCoV-2019_37'},
                        {'start': 11306, 'end': 11693, 'name': 'nCoV-2019_38'},
                        {'start': 11555, 'end': 11949, 'name': 'nCoV-2019_39'},
                        {'start': 11863, 'end': 12256, 'name': 'nCoV-2019_40'},
                        {'start': 12110, 'end': 12490, 'name': 'nCoV-2019_41'},
                        {'start': 12417, 'end': 12802, 'name': 'nCoV-2019_42'},
                        {'start': 12710, 'end': 13096, 'name': 'nCoV-2019_43'},
                        {'start': 13005, 'end': 13400, 'name': 'nCoV-2019_44'},
                        {'start': 13307, 'end': 13699, 'name': 'nCoV-2019_45_alt2'},
                        {'start': 13599, 'end': 13984, 'name': 'nCoV-2019_46'},
                        {'start': 13918, 'end': 14299, 'name': 'nCoV-2019_47'},
                        {'start': 14207, 'end': 14601, 'name': 'nCoV-2019_48'},
                        {'start': 14545, 'end': 14926, 'name': 'nCoV-2019_49'},
                        {'start': 14865, 'end': 15246, 'name': 'nCoV-2019_50'},
                        {'start': 15171, 'end': 15560, 'name': 'nCoV-2019_51'},
                        {'start': 15481, 'end': 15886, 'name': 'nCoV-2019_52'},
                        {'start': 15827, 'end': 16209, 'name': 'nCoV-2019_53'},
                        {'start': 16118, 'end': 16510, 'name': 'nCoV-2019_54'},
                        {'start': 16416, 'end': 16833, 'name': 'nCoV-2019_55'},
                        {'start': 16748, 'end': 17152, 'name': 'nCoV-2019_56'},
                        {'start': 17065, 'end': 17452, 'name': 'nCoV-2019_57'},
                        {'start': 17381, 'end': 17761, 'name': 'nCoV-2019_58'},
                        {'start': 17674, 'end': 18062, 'name': 'nCoV-2019_59'},
                        {'start': 17966, 'end': 18348, 'name': 'nCoV-2019_60'},
                        {'start': 18253, 'end': 18672, 'name': 'nCoV-2019_61'},
                        {'start': 18596, 'end': 18979, 'name': 'nCoV-2019_62'},
                        {'start': 18896, 'end': 19297, 'name': 'nCoV-2019_63'},
                        {'start': 19204, 'end': 19616, 'name': 'nCoV-2019_64'},
                        {'start': 19548, 'end': 19939, 'name': 'nCoV-2019_65'},
                        {'start': 19844, 'end': 20255, 'name': 'nCoV-2019_66'},
                        {'start': 20172, 'end': 20572, 'name': 'nCoV-2019_67'},
                        {'start': 20472, 'end': 20890, 'name': 'nCoV-2019_68'},
                        {'start': 20786, 'end': 21169, 'name': 'nCoV-2019_69'},
                        {'start': 21075, 'end': 21455, 'name': 'nCoV-2019_70'},
                        {'start': 21357, 'end': 21743, 'name': 'nCoV-2019_71'},
                        {'start': 21658, 'end': 22038, 'name': 'nCoV-2019_72'},
                        {'start': 21961, 'end': 22346, 'name': 'nCoV-2019_73'},
                        {'start': 22262, 'end': 22650, 'name': 'nCoV-2019_74'},
                        {'start': 22516, 'end': 22903, 'name': 'nCoV-2019_75'},
                        {'start': 22797, 'end': 23214, 'name': 'nCoV-2019_76'},
                        {'start': 23122, 'end': 23522, 'name': 'nCoV-2019_77'},
                        {'start': 23443, 'end': 23847, 'name': 'nCoV-2019_78'},
                        {'start': 23789, 'end': 24169, 'name': 'nCoV-2019_79'},
                        {'start': 24078, 'end': 24467, 'name': 'nCoV-2019_80'},
                        {'start': 24391, 'end': 24789, 'name': 'nCoV-2019_81'},
                        {'start': 24696, 'end': 25076, 'name': 'nCoV-2019_82'},
                        {'start': 24978, 'end': 25369, 'name': 'nCoV-2019_83'},
                        {'start': 25279, 'end': 25673, 'name': 'nCoV-2019_84'},
                        {'start': 25601, 'end': 25994, 'name': 'nCoV-2019_85'},
                        {'start': 25902, 'end': 26315, 'name': 'nCoV-2019_86'},
                        {'start': 26197, 'end': 26590, 'name': 'nCoV-2019_87'},
                        {'start': 26520, 'end': 26913, 'name': 'nCoV-2019_88'},
                        {'start': 26835, 'end': 27227, 'name': 'nCoV-2019_89'},
                        {'start': 27141, 'end': 27533, 'name': 'nCoV-2019_90'},
                        {'start': 27446, 'end': 27854, 'name': 'nCoV-2019_91'},
                        {'start': 27784, 'end': 28172, 'name': 'nCoV-2019_92'},
                        {'start': 28081, 'end': 28464, 'name': 'nCoV-2019_93'},
                        {'start': 28394, 'end': 28779, 'name': 'nCoV-2019_94'},
                        {'start': 28677, 'end': 29063, 'name': 'nCoV-2019_95'},
                        {'start': 28985, 'end': 29378, 'name': 'nCoV-2019_96'},
                        {'start': 29288, 'end': 29693, 'name': 'nCoV-2019_97'},
                        {'start': 29486, 'end': 29866, 'name': 'nCoV-2019_98'}];
                
                const svg = d3.select(selector)
                    .attr("width", svgwidth)
                    .attr("height", svgheight);
                
                const width = svgwidth;
                const height = svgheight;
                const title = 'Allele frequency: ' + sample;
                const xValue = d => d.POS;
                const xAxisLabel = 'Position';
                const yValue = d => d.FREQ;
                const yAxisLabel = 'Frequency';
                const y2AxisLabel = 'Coverage';
                const circleRadius = 4;
                const strokeWidth = 1.5;
                const legendRadiusSize = 5;
                const legendSpacing = 4;
                
                const margin = {top: 60, right: 190, bottom: 88, left: 80};
                const innerWidth = width - margin.left - margin.right;
                const innerHeight = height - margin.top - margin.bottom;

                const xScale = d3.scaleLinear()
                                .domain([0, d3.max(covdata, d => d.POS)])
                                .range([0, innerWidth])
                                .nice();
                
                const yScale = d3.scaleLinear()
                                .domain(d3.extent(data, yValue))
                                .range([innerHeight, 0])
                                .nice();
                
                const yScalecov = d3.scaleLinear()
                                .domain(d3.extent(covdata, d => d.DP))
                                .range([innerHeight, 0])
                                .nice();

                const yAxis = d3.axisLeft(yScale)
                            .tickSize(-innerWidth)
                            .ticks(10);
                
                const yAxiscov = d3.axisRight(yScalecov)
                            .tickSize(innerWidth)
                            .ticks(10);

                const xAxis = d3.axisBottom(xScale)
                            .tickSize(-innerHeight)
                            .tickPadding(20);

                //to set axis
                const g = svg.append('g')
                            .attr('transform', `translate(${margin.left}, ${margin.top})`);
                
                //Add labels and ticks to axes and remove some ticks
                const yAxisG = g.append('g')
                            .attr("class", "axisGrey")
                            .call(yAxis);
                
                /*
                .append("text")
                            .attr("transform", "rotate(-90)")
                            .attr("y", 6)
                            .attr("dy", ".71em")
                            .style("text-anchor", "end")
                            .text("Frequency (%))")
                */

                yAxisG.selectAll('.domain')
                        .remove();

                const xAxisG = g.append('g')
                            .attr("class", "axisGrey")
                            .call(xAxis)
                            .attr('transform', `translate(0, ${innerHeight})`);
                        
                xAxisG.select('.domain')
                    .remove();
                
                const yAxiscovG = g.append('g')
                            .attr("class", "axisNone")
                            .call(yAxiscov);
                
                yAxiscovG.select('.domain')
                    .remove();
                

                //Add x label
                xAxisG.append('text')
                    .attr('class', 'xAxisLabel')
                    .attr('y',60)
                    .attr('x', innerWidth / 2)
                    .attr('fill', 'black')
                    .text(xAxisLabel);

                //Add y label
                yAxisG.append('text')
                    .attr('class', 'yAxisLabel')
                    .attr('y', -50)
                    .attr('x', - innerHeight / 2)
                    .attr('fill', 'black')
                    .attr('transform', `rotate(-90)`)
                    .attr('text-anchor', 'middle')
                    .text(yAxisLabel);
                
                //Add y2 label
                yAxisG.append('text')
                    .attr('class', 'yAxisLabel')
                    .attr('y', innerWidth + 100)
                    .attr('x', - innerHeight / 2)
                    .attr('fill', 'black')
                    .attr('transform', `rotate(-90)`)
                    .attr('text-anchor', 'middle')
                    .text(y2AxisLabel);

                //Add title
                g.append('text')
                    .attr('class', 'title')
                    .attr('x', innerWidth / 2 )
                    .attr('y', -10)
                    .text(title);
                
                //Colors 
                //https://github.com/d3/d3-color/issues/11
                var lookup = {};
                var uniqueValues = [];

                data.forEach(d => {
                    if (!(d.ANN in lookup)) {
                        lookup[d.ANN] = 1;
                        uniqueValues.push(d.ANN);
                    };
                });
                
                var color = d3.scaleOrdinal()
                .domain(uniqueValues)
                .range(uniqueValues.map((group, i) => {
                    const t = i / uniqueValues.length;

                    return d3.hcl(t * 360, 50, 80);
                }));
            
                
                //Background rect primers
                g.selectAll('rect')
                .data(beddata)
                .enter()
                .append('rect')
                .attr('class', 'amplicon')
                .attr('x', d => xScale(d.start))
                .attr('y', d => yScale())
                .attr('width', d => xScale(d.end) - xScale(d.start))
                .attr('height', innerHeight)
                .attr('fill', '#fcd47e')
                .attr('fill-opacity', 0.05)
                .append('title')
                    .text(d => d.name + "\\n" + d.start + ":" + d.end);
                
                //Create line
                
                const line = d3.line()
                    .x(function(d) { return xScale(d.POS); })
                    .y(function(d) { return yScalecov(d.DP); });
                
                g.append("path")
                    .datum(covdata)
                    .attr("class", "line")
                    .attr("opacity", 0.8)
                    .attr("d", line);
                
                // Define the div for the frequency tooltip
                var div = d3.select("body").append("div")
                    .attr("class", "tooltip")
                    .style("opacity", 0);
                
                //Create scatter
                g.selectAll('circle')
                    .data(data)
                    .enter()
                    .append('circle')
                    .attr('fill-opacity', 0.8)
                    .attr('stroke-opacity', 1)
                    .attr('stroke', d => color(d.ANN))
                    .attr('fill', function(d) {return (yValue(d) > 50 ? "blue" : "darkseagreen"); })
                    .attr('stroke-width', strokeWidth)
                    .attr('cy', d => yScale(yValue(d)))
                    .attr('cx', d => xScale(xValue(d)))
                    .attr('r', circleRadius)
                    .on("mouseover", function(d) {
                        d3.select(this)
                        .transition()
                        .duration(50)
                        .attr('stroke-width', strokeWidth + 1);
                        div.transition()
                            .duration(100)
                            .style("opacity", .9);
                        div.html('POS: ' + d.POS + "<br/>" 
                                + 'FREQ: ' + d.FREQ + "<br/>" 
                                + d.REF + ">" + d.ALT + "<br/>"
                                + 'ANN: ' + d.ANN)
                            .style("left", (d3.event.pageX + 3) + "px")
                            .style("top", (d3.event.pageY - 30) + "px");
                        })
                    .on("mouseout", function(d) {
                        d3.select(this)
                        .transition()
                        .duration(500)
                        .attr('stroke-width',strokeWidth)

                        div.transition()
                            .duration(500)
                            .style("opacity", 0)
                    });

                //Add legend
                
                var legendG = svg.append("g")
                    .attr("class", "legend")
                    .attr('transform', `translate(${width - (margin.right / 2) - 30}, ${height - (margin.bottom * 2) + 30})`)
                    .attr("height", 100)
                    .attr("width", 100);
                
                var legend = legendG.selectAll('.legend')                     
                .data(color.domain())                                   
                .enter()                                                
                .append('g')                                            
                .attr('class', 'legend') 
                .attr('transform', function(d, i) {                 
                    var height = legendRadiusSize + (legendRadiusSize * 2);          
                    var offset =  height * color.domain().length / 2;     
                    var horz = -2 * legendRadiusSize;                       
                    var vert = i * height - offset;                       
                    return 'translate(' + horz + ',' + vert + ')';        
                });
                legend.append('circle')
                    .attr('id', d => d)
                    .attr('r', legendRadiusSize)
                    .attr('fill', d => color(d))
                    .style('stroke', d => color(d));
                
                legend.append('text')
                    .attr('x', legendRadiusSize + legendSpacing)
                    .attr('y', legendRadiusSize - legendSpacing)
                    .attr('dy', '0.32em')
                    .text(d => d);

                };

        DATAGROUP
                
        RENDERGROUP
                
        </script>
                    
        </body>
    </html>
    """



def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all files')
        parser.add_argument('-o', '--output', type=str, required=False, help='Output file to extract all results')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)

    if args.output == None:
        output_file = os.path.join(input_dir, 'vcf_report.html')
    else:
        output_file = args.output

    output_file = os.path.abspath(output_file)
    output_dir = ('/').join(output_file.split("/")[0:-1])

    check_create_dir(output_dir)

    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.date.today())
    right_now_full = "_".join(right_now.split(" "))

    log_filename = 'draw_vcf' + "_" + right_now_full + ".log"
    log_full_path = os.path.join(output_dir, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    #stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    
    #####################START PIPELINE################

    logger.info(args)
    
    vcf_files, cov_files = find_files(input_dir)

    if len(vcf_files) != len(cov_files):
        logger.error('Some missing files')
        sys.exit(1)

    svg_group = ""
    data_group = ""
    render_group = ""
    for vcf, cov in zip(vcf_files, cov_files):
        sample = vcf.split("/")[-1].split(".")[0]
        
        dfcov = import_cov_to_pandas(cov)
        dfvcf = import_VCF41_to_pandas(vcf)
        
        json_cov = dfcov.to_dict(orient='records')
        json_vcf = dfvcf.to_dict(orient='records')
    
    
        svg_group = svg_group + '<svg id="S' + sample  + '\"></svg>' + "\n"
        
        data_group = data_group + "const dataVCF_" + sample + " = " + str(json_vcf) + ";\n" + \
                    "const dataCOV_" + sample + " = " + str(json_cov) + ";\n"
        
        render_group = render_group + 'render("#S' + sample + '", dataVCF_' + sample + \
                    ', dataCOV_' + sample + ', "' + sample + '")' + "\n"



    template_file = d3_template.replace('SVGGROUP', svg_group)
    template_file = template_file.replace('DATAGROUP', data_group)
    template_file = template_file.replace('RENDERGROUP', render_group)
    with open(output_file,'w+') as fout:
        fout.write(template_file)
    


    logger.info("DONE: file can be found in " + output_file)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise