$(document).ready( function () {
    $('#probstable').DataTable({
        dom: '<"dom_wrapper fh-fixedHeader"Bf>tip',
        buttons: ['copy', 'excel', 'pdf' ],
        fixedHeader: true,
        order: [[ 3, "desc" ]],
        columnDefs: [
            {
                targets: [0],
                className: "dt-left",
                render: function (data, type, row, meta) {
                    if(type === 'display') {
                        data = '<a target="_blank" href="https://www.ncbi.nlm.nih.gov/gene/' + encodeURIComponent(data) + '">' + data + '</a>';
                    }
                    return data;
                }

            },
            {
                targets: [1, 2],
                className: "dt-left"
            },
            {
                targets: [ 3 ],
                className: "dt-right",
                defaultOrder: true,
                sortOrder: 'desc',
                render: function (data, type, full) {
                    return parseFloat(data).toFixed(2);
                }
            },
            {
                targets: [ 4 ],
                className: "dt-center"
            },
            {
                targets: [ 5 ],
                className: "dt-center"
            }  ,
            {
                targets: [ 6 ],
                className: "dt-right"
            }
        ],
        pagingType: "full_numbers",
        initComplete: function(){
            $("#probstable").show();
        }
    });

    $('#gotable').DataTable({
        responsive: {
            details: {
                display: $.fn.dataTable.Responsive.display.modal( {
                    header: function ( row ) {
                        var data = row.data();
                        return 'Details for '+data[0]+' '+data[1];
                    }
                } ),
                renderer: $.fn.dataTable.Responsive.renderer.tableAll({
                    tableClass: 'table'
                })
            }
        },
        dom: '<"dom_wrapper fh-fixedHeader"Bf>tip',
        buttons: ['copy', 'excel', 'pdf' ],
        fixedHeader: true,
        order: [[ 2, "desc" ]],
        columnDefs: [
            {
                responsivePriority: 1,
                targets: [0],
                className: "dt-left all",
                render: function (data, type, row, meta) {
                    if(type === 'display') {
                        data = '<a target="_blank" href="http://amigo.geneontology.org/amigo/term/' + encodeURIComponent(data) + '">' + data + '</a>';
                    }
                    return data;
                }
            },
            {
                targets: [1],
                className: "dt-left all"
            },
            {
                targets: [2],
                className: "dt-right all",
                defaultOrder: true,
                sortOrder: 'desc',
                render: function (data, type, full) {
                    return parseFloat(data).toFixed(2);
                }
            },
            {
                targets: [3],
                className: "dt-right all"

            }
        ],
        pagingType: "full_numbers",
        initComplete: function(){
            $("#gotable").show();
        }
    });

    $('#distable').DataTable({
        dom: '<"dom_wrapper fh-fixedHeader"Bf>tip',
        buttons: ['copy', 'excel', 'pdf' ],
        fixedHeader: true,
        order: [[ 2, "desc" ]],
        columnDefs: [
            {
                targets: [0],
                className: "dt-left",
                render: function (data, type, row, meta) {
                    if(type === 'display') {
                        data = '<a target="_blank" href="https://disease-ontology.org/?id=' + encodeURIComponent(data) + '">' + data + '</a>';
                    }
                    return data;
                }
            },
            {
                targets: [1],
                className: "dt-left"
            },
            {
                targets: [2],
                className: "dt-right",
                defaultOrder: true,
                sortOrder: 'desc',
                render: function (data, type, full) {
                    return parseFloat(data).toFixed(2);
                }
            },
            {
                targets: [3],
                className: "dt-right all"

            }
        ],
        pagingType: "full_numbers",
        initComplete: function(){
            $("#distable").show();
        }
    });


    $('#validatetable').DataTable({
        dom: '<"dom_wrapper fh-fixedHeader"Bf>tip',
        columnDefs: [
            {
                targets: [0, 1],
                className: "dt-left"
            },
            {
                targets: [2, 3, 4, 5],
                className: "dt-center"
            }
        ],
        pagingType: "full_numbers",
        buttons: ['copy', 'excel', 'pdf' ],
        initComplete: function(){
            $("#validatetable").show();
        }
    });

    $('#validateresults').DataTable({
        dom: '<"dom_wrapper fh-fixedHeader"Bf>tip',
        columnDefs: [
            {
                targets: [0, 1],
                className: "dt-left"
            },
            {
                targets: [2],
                className: "dt-center"
            }
        ],
        pagingType: "full_numbers",
        buttons: ['copy', 'excel', 'pdf' ],
        initComplete: function(){
            $("#validateresults").show();
        }
    });

} );