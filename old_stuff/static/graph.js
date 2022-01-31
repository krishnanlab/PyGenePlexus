//Create SVG element
var svg = d3.select("svg"),
    w = +svg.node().getBoundingClientRect().width,
    h = +svg.node().getBoundingClientRect().height

//////// g holds zoom, nodes and links
var g = svg.append('g')
    .attr("class", "everything");

// define scale variable and range for node sizing
var nodescale = d3.scaleLinear().domain([0, .3]).range([1, 10])


var items = [
    {
      label: function (d) {
        return d.id;
      },
      items: [
        {
          label: "Entrez Lookup",
          onClick: function (d) {
          console.log(d.id)
            window.open('https://www.ncbi.nlm.nih.gov/gene/' + d.id, '_blank');

          }
        }
      ]
    }
  ];

/*
//manipulate stuff
d3.selectAll(".resize").append("circle")
    .attr("cx", 0)
    .attr("cy", 0)
    .attr("r", radius)
    .attr("fill", function(d){return myColor(d.Class) })
    .classed("handle", true)


d3.select(".domain")
  .select(function () { return this.parentNode.appendChild(this.cloneNode(true)) })
  .classed("halo", true)
*/

// Setup the tool tip.
var tool_tip = d3.tip()
  .attr('class', 'd3-tip')
  .offset([-10, 0])
  .html(function(d) {
      return "<strong>Symbol: " + d.Symbol + "</strong><br>" +
      "<strong>Entrez: " + d.id + "</strong><br>" +
      "<strong>Probability: " + d.Probability + "</strong><br>" +
      "<strong>Class-Label: " + d.Class + "</strong><br>";
  });

svg.call(tool_tip);

// Color the nodes by Class
var data = ['P','U','N']
var myColor = d3.scaleOrdinal().domain(data)
    .range(["#00ccbc","#2c7bb6","#CC333F"]);


//////////// FORCE SIMULATION ////////////
var simulation = d3.forceSimulation();

// set up the simulation and event to update locations after each tick
function initializeSimulation() {
  simulation.nodes(dataset.nodes);
  initializeForces();
  simulation.on("tick", ticked);
}

// values for all forces
forceProperties = {
    center: {
        x: 0,
        y: 0
    },
    charge: {
        enabled: true,
        strength: -1500,
        distanceMin: 1,
        distanceMax: 150
    },
    collide: {
        enabled: true,
        strength: .7,
        iterations: 1,
        radius: 5
    },
    link: {
        enabled: true,
        distance: 100,
        iterations: 1
    }
}

// add forces to the simulation
function initializeForces() {
    // add forces and associate each with a name
    simulation
        .force("link", d3.forceLink())
        .force("charge", d3.forceManyBody())
        .force("collide", d3.forceCollide())
        .force("center", d3.forceCenter())
    // apply properties to each of the forces
    updateForces();
}

// apply new force properties
function updateForces() {
    simulation.force("charge")
        .strength(forceProperties.charge.strength * forceProperties.charge.enabled)
        .distanceMin(forceProperties.charge.distanceMin)
        .distanceMax(forceProperties.charge.distanceMax);
    simulation.force("collide")
        .strength(forceProperties.collide.strength * forceProperties.collide.enabled)
        .radius(forceProperties.collide.radius)
        .iterations(forceProperties.collide.iterations);
    simulation.force("link")
         .id(function(d) {return d.id})
        .distance(forceProperties.link.distance)
        .iterations(forceProperties.link.iterations)
        .links(forceProperties.link.enabled ? dataset.links : []);

    // updates ignored until this is run
    // restarts the simulation (important if simulation has already slowed down)
    simulation.alpha(1).restart();
}

//////////// DISPLAY ////////////

// generate the svg objects and force simulation
function initializeDisplay() {

    var menu = contextMenu().items('Entrez Lookup', 'Symbol Lookup');

    var toggle = 0;
    // functions for highlighting neighboring links
    /*
    const linkedByIndex = {};
    dataset.links.forEach(d => {
        linkedByIndex[`${d.source.index},${d.target.index}`] = 1;
    });

    function isConnected(a, b) {
        return linkedByIndex[`${a.index},${b.index}`] || linkedByIndex[`${b.index},${a.index}`] || a.index === b.index;
    }

    function fade(opacity) {
        return d => {
          node.style('stroke-opacity', function (o) {
            const thisOpacity = isConnected(d, o) ? 1 : opacity;
            this.setAttribute('fill-opacity', thisOpacity);
            return thisOpacity;
          });

          link.style('stroke-opacity', o => (o.source === d || o.target === d ? 1 : opacity));

        };
    }

    // New fade functions
    var linkedByIndex = {};
    dataset.links.forEach(function(d) {
      linkedByIndex[d.source.index + "," + d.target.index] = 1;
    });

    function isConnected(a, b) {
        return linkedByIndex[a.index + "," + b.index] || linkedByIndex[b.index + "," + a.index] || a.index == b.index;
    }

    function fade(opacity) {
        return function(d) {
            node.style("stroke-opacity", function(o) {
                thisOpacity = isConnected(d, o) ? 1 : opacity;
                this.setAttribute('fill-opacity', thisOpacity);
                return thisOpacity;
            });

            link.style("stroke-opacity", opacity).style("stroke-opacity", function(o) {
                return o.source === d || o.target === d ? 1 : opacity;
            });
        };
    }
    */
    // set the data and properties of link lines
    link = g.append("g")
        .attr("class", "links")
        .selectAll("line")
        .data(dataset.links)
        .enter().append("line")
        .style("stroke", "#ADA9A8")
        .style("stroke-width", function(d) { return (d.weight); })
        //.on('mouseout.fade', fade(1));


    // set the data and properties of node circles
    node = g.append("g")
        .attr("class", "nodes")
        .selectAll("circle")
        .data(dataset.nodes)
        .enter()
        // add note g element for each node here.
        .append("g")
        // position the g element like the circle element use to be.


    node.append("circle")
        .attr("r", function(d){return nodescale(d.Probability)})
        //.attr("r", function(d){ return Math.exp(d.Probability)*10})
        .attr("fill", function(d){return myColor(d.Class) })
        .on('mouseover', tool_tip.show)
        //.on('mouseover.fade', fade(0.1))
        .on('mouseout', tool_tip.hide)
  	    //.on('mouseout.fade', fade(1))
        .on('contextmenu', d3.contextmenu(items));

    node.append("text")
        .attr("text-anchor", "middle")
        .text(function(d) { return d.Symbol; });

    //add drag capabilities
    var drag_handler = d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended);

    drag_handler(node);

    //add zoom capabilities
    var zoom_handler = d3.zoom()
        .on("zoom", zoom_actions);

    zoom_handler(svg);

  // visualize the graph
  updateDisplay();
}

// update the display based on the forces (but not positions)
function updateDisplay() {
    node
        //.attr("stroke", forceProperties.charge.strength > 0 ? "blue" : "red")
        .attr("stroke-width", forceProperties.charge.enabled==false ? 0 : Math.abs(forceProperties.charge.strength)/15);

    link
        //.attr("stroke-width", forceProperties.link.enabled ? 1 : (d.weight))
        .attr("opacity", forceProperties.link.enabled ? 1 : 0);
}

// update the display positions after each simulation tick
function ticked() {
    link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node
        //.attr("cx", function(d) { return d.x; })
        //.attr("cy", function(d) { return d.y; })
        .attr("transform", function(d) {
          return "translate(" + d.x + "," + d.y + ")";
        })

    d3.select('#alpha_value').style('flex-basis', (simulation.alpha()*100) + '%');
}


//////////// UI EVENTS ////////////

//Drag functions
//d is the node
function dragstarted(d) {
  if (!d3.event.active) simulation.alphaTarget(0.3).restart();
  d.fx = d.x;
  d.fy = d.y;
}
//make sure you can't drag the circle outside the box
function dragged(d) {
  d.fx = d3.event.x;
  d.fy = d3.event.y;
}

function dragended(d) {
  if (!d3.event.active) simulation.alphaTarget(0.0001);
  d.fx = null;
  d.fy = null;
}

//Zoom functions
function zoom_actions(){
    g.attr("transform", d3.event.transform)
}


//Context Menu functions
function contextMenu() {
    var height,
        width,
        margin = 0.1, // fraction of width
        items = [],
        style = {
            'rect': {
                'mouseout': {
                    'fill': 'rgb(244,244,244)',
                    'stroke': 'white',
                    'stroke-width': '1px'
                },
                'mouseover': {
                    'fill': 'rgb(200,200,200)'
                }
            },
            'text': {
                'fill': 'steelblue',
                'font-size': '13'
            }
        };

    function menu(x, y) {
        d3.select('.context-menu').remove();

        // Draw the menu
        d3.select('g')
            .append('g').attr('class', 'context-menu')
            .selectAll('tmp')
            .data(items).enter()
            .append('g').attr('class', 'menu-entry')
            .style({'cursor': 'pointer'})
            .on('mouseover', function(){
                d3.select(this).select('rect').style(style.rect.mouseover) })
            .on('mouseout', function(){
                d3.select(this).select('rect').style(style.rect.mouseout) });

        d3.selectAll('.menu-entry')
            .append('rect')
            .attr('x', x)
            .attr('y', function(d, i){ return y + (i * height); })
            .attr('width', width)
            .attr('height', height)
            .style(style.rect.mouseout);

        d3.selectAll('.menu-entry')
            .append('text')
            .text(function(d){ return d; })
            .attr('x', x)
            .attr('y', function(d, i){ return y + (i * height); })
            .attr('dy', height - margin / 2)
            .attr('dx', margin)
            .style(style.text);

        // Other interactions
        d3.select('body')
            .on('click', function() {
                d3.select('.context-menu').remove();
            });

    }

    menu.items = function(e) {
        if (!arguments.length) return items;
        for (i in arguments) items.push(arguments[i]);
        return menu;
    }

    return menu;
}



// update size-related forces
d3.select(window).on("resize", function(){
    width = +svg.node().getBoundingClientRect().width;
    height = +svg.node().getBoundingClientRect().height;
    updateForces();
});

// convenience function to update everything (run after UI input)
function updateAll() {
    updateForces();
    updateDisplay();
}

initializeDisplay();
initializeSimulation();