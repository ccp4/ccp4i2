import { Card, CardHeader, Select } from '@mui/material'
import { Fragment, useEffect, useState } from 'react'
import { ComposedChart, CartesianGrid, XAxis, YAxis, Legend, Tooltip, Scatter, Bar, Label, ReferenceArea } from 'recharts'
import $ from 'jquery'

export const reportXmlToCharts = (reportXML) => {
    const decodeHTMLEntities = (str) => {
        if (str && typeof str === 'string') {
            str = str.replace(/&nbsp;/g, ' ');
        }
        return str;
    }
    var newCharts = []
    var newTables = $(reportXML).find('ccp4_data').toArray()
    if (newTables.length === 0) {
        newTables = $(reportXML).find('ccp4\\:ccp4_data').toArray()
    }
    newTables.forEach(item => {
        var newChart = { plots: [], datas: {} }
        var headersNode = $(item).find('headers')
        var separator = headersNode.attr('separator')
        if (separator == null || separator === ' ') separator = /\s+|\t+/
        newChart.headers = headersNode.text().trim().split(separator).map(item => decodeHTMLEntities(item))
        $(item).find('data').toArray().forEach((data, iData) => {
            var dataId = $(data).attr('id')
            if (dataId == null) {
                dataId = `${iData}`
            }
            newChart.datas[dataId] = $(data).text().trim().split(/\n|\r/)
                .map((rowString, iRowString) => {
                    var row = rowString.trim().split(/\s+|\t+/)
                    var rowDict = { 0: iRowString }
                    row.forEach((item, iItem) => {
                        rowDict[iItem + 1] = parseFloat(row[iItem])
                        if (Number.isNaN(rowDict[iItem + 1])) {
                            rowDict[iItem + 1] = row[iItem]
                        }
                    })
                    return rowDict
                })
            //console.log('Interim keys', Object.keys(newChart.datas))
        })
        //console.log('Initially keys', Object.keys(newChart.datas))
        var chartPlots = $(item).find('plot').toArray()
        chartPlots.forEach(plot => {
            var newPlot = {}
            newPlot.plottype = $(plot).find('plottype').text().trim()
            newPlot.title = $(plot).find('title').text().trim()
            newPlot.xLabel = decodeHTMLEntities($(plot).find('xlabel').text().trim())
            newPlot.yLabel = decodeHTMLEntities($(plot).find('ylabel').text().trim())
            newPlot.ryLabel = decodeHTMLEntities($(plot).find('rylabel').text().trim())
            newPlot.plotLines = []
            newPlot.barCharts = []
            newPlot.yAxes = {}//{ left: { yAxisId: 'left', label: '' }, right: { yAxisId: 'right', label: '', orientation: 'right' } }
            newPlot.xAxes = {}
            newPlot.leftYAxes = {}
            newPlot.rightYAxes = {}
            newPlot.xintegral = $(plot).find('xintegral').text()
            newPlot.yintegral = $(plot).find('yintegral').text()
            newPlot.xscale = $(plot).find('xscale').text()
            newPlot.polygons = []
            newPlot.circles = []
            newPlot.lines = []
            newPlot.fixaspectratio = false

            newPlot.xrange = $(plot).find('xrange')
            newPlot.xrange = { min: parseFloat($(newPlot.xrange).attr('min')), max: parseFloat($(newPlot.xrange).attr('max')) }
            if (newPlot.xrange.min == null || Number.isNaN(newPlot.xrange.min)) newPlot.xrange.min = 'auto'
            if (newPlot.xrange.max == null || Number.isNaN(newPlot.xrange.max)) newPlot.xrange.max = 'auto'
            if ($(plot).find("fixaspectratio").text() === "true") { newPlot.fixaspectratio = true }
            //console.log('xrange', newPlot.xrange)
            newPlot.yrange = $(plot).find('yrange')
            newPlot.yrange = { min: parseFloat($(newPlot.yrange).attr('min')), max: parseFloat($(newPlot.yrange).attr('max')) }
            if (newPlot.yrange.min == null || Number.isNaN(newPlot.yrange.min)) newPlot.yrange.min = 'auto'
            if (newPlot.yrange.max == null || Number.isNaN(newPlot.yrange.max)) newPlot.yrange.max = 'auto'

            //console.log('yrange', newPlot.yrange)
            var plotLines = $(plot).find('plotline').toArray()
            var barCharts = $(plot).find('barchart').toArray()
            let allTypes = plotLines.concat(barCharts)

            //console.log('newPlot', newPlot)
            //Handle polygons
            newPlot.polygons = $(plot).find('polygon').toArray().map((polygon, iPolygon) => {
                var newPolygon = {
                    points: $(polygon).text().trim().split(/\s+|\t+/),
                    fill: $(polygon).attr("fillcolour"),
                    stroke: $(polygon).attr("linecolour"),
                    fillOpacity: $(polygon).attr("alpha"),
                }
                return newPolygon
            })
            //Handle circles
            newPlot.circles = $(plot).find('circle').toArray().map((circle, iCircle) => {
                var newCircle = {
                    cx: parseFloat($(circle).attr("xpos")),
                    cy: parseFloat($(circle).attr("ypos")),
                    r: parseFloat($(circle).attr("radius")),
                    stroke: $(circle).attr("linecolour"),
                    fillOpacity: parseFloat($(circle).attr("alpha")),
                }
                if (newCircle.fillOpacity == null || Number.isNaN(newCircle.fillOpacity)) {
                    newCircle.fillOpacity = 0.
                }
                return newCircle
            })
            //Handle lines
            //<line x1="  -4.0000" x2="   4.0000" y1="  -4.0000" y2="   4.0000" linestyle="-" linecolour="black"/>
            newPlot.lines = $(plot).find('line').toArray().map((line, iLine) => {
                var newLine = {
                    x1: parseFloat($(line).attr("x1")),
                    x2: parseFloat($(line).attr("x2")),
                    y1: parseFloat($(line).attr("y1")),
                    y2: parseFloat($(line).attr("y2")),
                    linestyle: $(line).attr("linestyle"),
                    stroke: $(line).attr("linecolour"),
                }
                return newLine
            })

            allTypes.forEach((plotLine, iPlotLine) => {
                var newPlotLine = {
                    type: plotLine.nodeName,
                    xcol: parseInt($(plotLine).attr('xcol')),
                    ycol: parseInt($(plotLine).attr('ycol')),
                    dataid: $(plotLine).attr('ycol'),
                    color: $(plotLine).attr('colour'),
                    linestyle: $(plotLine).find('linestyle').text(),
                    label: $(plotLine).find('label').text(),
                    symbolsize: parseInt($(plotLine).find("symbolsize").text())
                }
                //console.log(newPlotLine.symbolsize)
                if (newPlotLine.symbolsize == null || Number.isNaN(newPlotLine.symbolsize)) {
                    newPlotLine.symbolsize = 1
                }
                if (newPlotLine.label == null || (typeof newPlotLine.color === "label" && newPlotLine.label.length == 0)) {
                    newPlotLine.label = false
                }
                if (plotLine.nodeName === 'barchart') {
                    newPlotLine.xcol = parseInt($(plotLine).attr('col'))
                    newPlotLine.ycol = parseInt($(plotLine).attr('tcol'))
                }
                //console.log('plotLine had dataid', newPlotLine.dataid, Object.keys(newChart.datas))
                if (newPlotLine.dataid == null || !Object.keys(newChart.datas).includes(newPlotLine.dataid)) {
                    newPlotLine.dataid = Object.keys(newChart.datas)[0]
                }
                //console.log('plotLine has dataid', newPlotLine.dataid)

                if (newPlotLine.color == null || (typeof newPlotLine.color === "string" && newPlotLine.color.length == 0)) {
                    newPlotLine.color = $(plotLine).find("colour").text()
                    if (newPlotLine.color == null || (typeof newPlotLine.color === "string" && newPlotLine.color.length == 0)) {
                        newPlotLine.color = ["red", "green", "blue", "cyan", "magenta", "yellow"][iPlotLine % 6]
                    }
                }

                var xAxisId = newPlotLine.xcol//$(plotLine).attr('topaxis') === "true" ? "top" : "bottom"                    
                if (!Object.keys(newPlot.xAxes).includes(newPlotLine.xcol)) {
                    newPlot.xAxes[xAxisId] = {
                        plotLine: newPlotLine,
                        xAxisId: xAxisId,
                        dataKey: newPlotLine.xcol,
                        position: "insideBottom",
                        domain: [newPlot.xrange.min, newPlot.xrange.max]
                    }
                }
                if (newPlot.xscale === "oneoversqrt") {
                    newPlot.xAxes[xAxisId].tickFormatter = (value) => { return Math.pow(value, -0.5).toFixed(2) }
                }
                newPlotLine.xAxisId = xAxisId

                var yAxisId = newPlotLine.ycol//$(plotLine).attr('rightaxis') === "true" ? "right" : "left"
                if (!Object.keys(newPlot.yAxes).includes(yAxisId)) {
                    newPlot.yAxes[yAxisId] = {
                        plotLine: newPlotLine,
                        yAxisId: yAxisId,
                        dataKey: newPlotLine.ycol,
                        position: $(plotLine).attr('rightaxis') === "true" ? "insideRight" : "insideLeft",
                        orientation: $(plotLine).attr('rightaxis') === "true" ? "right" : "left",
                        angle: -90.,
                        domain: [newPlot.yrange.min, newPlot.yrange.max]
                    }
                }
                newPlotLine.yAxisId = yAxisId

                newPlot.plotLines.push(newPlotLine)
            })
            /*
            *Now some fun and games to make all axes span the same domain
            */
            Object.keys(newPlot.yAxes).forEach((axisKey, iAxisKey) => {
                var axis = newPlot.yAxes[axisKey]
                if (axis.orientation == "right") { newPlot.rightYAxes[axisKey] = axis }
                else { newPlot.leftYAxes[axisKey] = axis }
            })
            var axisSets = [newPlot.xAxes, newPlot.leftYAxes, newPlot.rightYAxes]

            axisSets.forEach((axisSet, iAxisSet) => {
                var globalMin = 1e30
                var globalMax = -1e30
                Object.keys(axisSet).forEach((axisKey, iAxis) => {
                    var axis = axisSet[axisKey]
                    var relevantData = newChart.datas[axis.plotLine.dataid].map(item => parseFloat(item[axis.dataKey]))
                    var filteredData = relevantData.filter(entry => !Number.isNaN(parseFloat(entry)))
                    var axisMin = axis.domain[0]
                    if (axisMin === 'auto') {
                        axisMin = Math.min(1e30, Math.min(...filteredData))
                    }
                    var axisMax = axis.domain[1]
                    if (axisMax === 'auto') {
                        axisMax = Math.max(-1e30, Math.max(...filteredData))
                    }
                    globalMin = globalMin < axisMin ? globalMin : axisMin
                    globalMax = globalMax > axisMax ? globalMax : axisMax
                })
                Object.keys(axisSet).forEach((axisKey, iAxis) => {
                    var axis = axisSet[axisKey]
                    axis.domain = [globalMin != null ? globalMin : "auto", globalMax != null ? globalMax : "auto"]
                })

            })


            newChart.plots.push(newPlot)


        })
        newCharts.push(newChart)
    })
    return newCharts
}

export const CCP4i2ReportChart = (props) => {
    const [selectedPlotOption, setSelectedPlotOption] = useState(0)
    const [newCharts, setNewChart] = useState([])
    const [plotOptions, setPlotOptions] = useState([])
    const [width, setWidth] = useState(0)
    const [height, setHeight] = useState(0)
    const [calculatedChart, setCalculatedChart] = useState([])

    useEffect(() => {
        if (props.chart) {
            console.log('chart', props.chart)
            setWidth(props.width)
            setHeight(props.height)
            var newOptions = props.chart.plots.map((plot, iItem) => { return { value: iItem, label: plot.title } })
            setPlotOptions(newOptions)

            var newCalculatedChart = []
            props.chart?.plots?.map((plot, iPlot) => {
                var newCalculatedPlot = {}
                var newScatters = []
                var newBarCharts = []

                plot.plotLines.map((plotLine, iPlotLine) => {
                    if (plotLine.type === "plotline") {
                        let params = {
                            key: `Scatter_${iPlotLine}`,
                            name: plotLine.label ? plotLine.label : props.chart.headers[parseInt(plotLine.ycol) - 1],
                            line: plotLine.linestyle === "." ? false : { stroke: plotLine.color, strokeWidth: 2 },
                            fill: plotLine.color,
                            xAxisId: plotLine.xAxisId,
                            yAxisId: plotLine.yAxisId,
                            dataKey: plotLine.ycol,
                            isAnimationActive: false,
                            shape: <MyCircle r={plotLine.symbolsize} />,
                            data: props.chart.datas[plotLine.dataid].filter(item => !(item[plotLine.xcol] === '-' || item[plotLine.ycol] === '-')),
                        }
                        if (plotLine.linestyle === ".") { delete params.dataKey }
                        else { delete params.data }
                        newScatters.push(<Scatter {...params} />)
                    }
                    else {
                        let params = {
                            key: `Bar_${iPlotLine}`,
                            yAxisId: plotLine.yAxisId,
                            xAxisId: plotLine.xAxisId,
                            name: props.chart.headers[plotLine.xcol],
                            dataKey: plotLine.ycol,
                            isAnimationActive: false,
                            fill: plotLine.color,
                            fillOpacity: 0.5,
                            stroke: plotLine.color,
                            strokeWidth: 2
                        }
                        newBarCharts.push(<Bar {...params} />)
                    }
                    newCalculatedPlot.scatters = newScatters
                    newCalculatedPlot.barCharts = newBarCharts

                    newCalculatedPlot.xAxes = Object.keys(plot.xAxes).map((axisId, iAxis) => {
                        var axis = plot.xAxes[axisId]
                        return <XAxis label={<Label position={axis.position}>{plot.xLabel}</Label>} key={`XAxis_${axisId}`}
                            {...axis}
                            type="number"
                            hide={iAxis != 0} />
                    })

                    newCalculatedPlot.leftYAxes = Object.keys(plot.leftYAxes).map((axisId, iAxis) => {
                        var axis = plot.yAxes[axisId]
                        delete axis.angle
                        return <YAxis label={<Label position={axis.position} angle={-90}>{axis.orientation === "left" ? plot.yLabel : plot.ryLabel}</Label>}
                            key={`LeftYAxis_${axisId}`} {...axis} type="number" hide={iAxis != 0}
                        />
                    })
                    newCalculatedPlot.rightYAxes = Object.keys(plot.rightYAxes).map((axisId, iAxis) => {
                        var axis = plot.yAxes[axisId]
                        delete axis.angle
                        return <YAxis label={<Label position={axis.position} angle={-90}>{axis.orientation === "left" ? plot.yLabel : plot.ryLabel}</Label>}
                            key={`RightYAxis_${axisId}`} {...axis} type="number" hide={iAxis != 0}
                        />
                    })
                    /*
                                    newCalculatedPlot.polygons = plot.polygons.map((polygon, iPolygon) => {
                                        return <polygon key={iPolygon} points={transformPoints(polygon.points, plot).join(",")} fill={polygon.fill} stroke={polygon.stroke} fill-opacity={polygon.fillOpacity} />
                                    })
                    */
                    newCalculatedPlot.polygons = plot.polygons.map((polygon, iPolygon) => {
                        var xMin = Math.min(...polygon.points.filter((item, iItem) => { return iItem % 2 == 0 }))
                        var xMax = Math.max(...polygon.points.filter((item, iItem) => { return iItem % 2 == 0 }))
                        var yMin = Math.min(...polygon.points.filter((item, iItem) => { return iItem % 2 == 1 }))
                        var yMax = Math.max(...polygon.points.filter((item, iItem) => { return iItem % 2 == 1 }))
                        var result = <ReferenceArea x1={xMin} y1={yMin} x2={xMax} y2={yMax}
                            stroke={polygon.stroke}
                            fill={polygon.fill}
                            fillOpacity={polygon.fillOpacity}
                            xAxisId={plot.xAxes[Object.keys(plot.xAxes)[0]].xAxisId}
                            yAxisId={plot.leftYAxes[Object.keys(plot.leftYAxes)[0]].yAxisId}
                            ifOverflow="extendDomain"
                            shape={<ReferenceRectangle />}
                        />
                        return result
                    })

                    newCalculatedPlot.circles = plot.circles.map((circle, iCircle) => {
                        var xMin = parseFloat(circle.cx) - parseFloat(circle.r)
                        var xMax = parseFloat(circle.cx) + parseFloat(circle.r)
                        var yMin = parseFloat(circle.cy) - parseFloat(circle.r)
                        var yMax = parseFloat(circle.cy) + parseFloat(circle.r)
                        var result = <ReferenceArea x1={xMin} y1={yMin} x2={xMax} y2={yMax}
                            stroke={circle.stroke}
                            xAxisId={plot.xAxes[Object.keys(plot.xAxes)[0]].xAxisId}
                            yAxisId={plot.leftYAxes[Object.keys(plot.leftYAxes)[0]].yAxisId}
                            ifOverflow="extendDomain"
                            shape={<ReferenceCircle />}
                        />
                        return result
                    })

                    newCalculatedPlot.lines = plot.lines.map((line, iLine) => {
                        var result = <ReferenceArea x1={line.x1} y1={line.y1} x2={line.x2} y2={line.y2}
                            stroke={line.stroke}
                            xAxisId={plot.xAxes[Object.keys(plot.xAxes)[0]].xAxisId}
                            yAxisId={plot.leftYAxes[Object.keys(plot.leftYAxes)[0]].yAxisId}
                            ifOverflow="extendDomain"
                            shape={<ReferenceLine />}
                        />
                        return result
                    })
                })
                newCalculatedChart.push(newCalculatedPlot)
            })
            setCalculatedChart(newCalculatedChart)
        }
    }, [props.chart])

    const transformPoints = (points, plot) => {
        var margins = [70 + Object.keys(plot.rightYAxes).length > 0 ? 62 : 0, 62]
        var offset = [65, 5]
        var transformedPoints = points.map((coord, iCoord) => {
            if (iCoord % 2 === 0) {
                var domain = plot.xAxes[Object.keys(plot.xAxes)[0]].domain
                return offset[0] + (coord - domain[0]) * ((props.width - margins[0]) / (domain[1] - domain[0]))
            }
            else {
                var domain = plot.yAxes[Object.keys(plot.yAxes)[0]].domain
                return offset[1] + (props.height - margins[1]) - ((coord - domain[0]) * ((props.height - margins[1]) / (domain[1] - domain[0])))
            }
        })
        return transformedPoints
    }

    const MyCircle = (arg) => {
        return <circle cx={arg.cx} cy={arg.cy} r={arg.r} stroke={arg.fill} fill="none" />
    }

    const ReferenceCircle = props => {
        return (
            <circle
                stroke={props.stroke}
                fill={props.fill}
                fillOpacity={0.2}
                cx={props.x + props.width / 2}
                cy={props.y + props.height / 2}
                r={props.width / 2}
            />
        );
    };

    const ReferenceLine = props => {
        return (
            <line
                stroke={props.stroke}
                x1={props.x} x2={props.x + props.width} y1={props.y + props.height} y2={props.y}
            />
        );
    };

    const ReferenceRectangle = props => {
        return (
            <polygon
                stroke={props.stroke}
                fill={props.fill}
                fillOpacity={props.fillOpacity}
                points={[
                    props.x, props.y,
                    props.x + props.width, props.y,
                    props.x + props.width, props.y + props.height,
                    props.x, props.y + props.height,
                    props.x, props.y].join(',')}
            />
        );
    };

    return <Card action={props.extra}>
        <CardHeader>
            {plotOptions.length > 1 ?
                <Fragment>Select plot: <Select value={selectedPlotOption} options={plotOptions}
                    onSelect={(value) => {
                        setSelectedPlotOption(value)
                        if (plotOptions[selectedPlotOption].fixaspectratio) { setHeight(props.width) }
                        else { setHeight(props.height) }
                    }} />
                </Fragment>
                : <h1>Plot title: {plotOptions.length > 0 && plotOptions[0].label}</h1>
            }
        </CardHeader>

        {calculatedChart.map((plot, iPlot) => {
            return iPlot === selectedPlotOption &&
                <ComposedChart width={width} height={height}
                    data={props.chart.datas[Object.keys(props.chart.datas)[0]]} key={`Plot_${iPlot}`}
                    margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" />
                    <Tooltip />
                    <Legend />
                    {plot.xAxes.map(item => item)}
                    {plot.leftYAxes.map(item => item)}
                    {plot.rightYAxes.map(item => item)}
                    {plot.scatters.map(item => item)}
                    {plot.barCharts.map(item => item)}
                    {plot.polygons.map(item => item)}
                    {plot.circles.map(item => item)}
                    {plot.lines.map(item => item)}
                </ComposedChart>
        })
        }
    </Card>
}
CCP4i2ReportChart.defaultProps = { width: 300, height: 300 }