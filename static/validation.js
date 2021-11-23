$('[data-toggle="popover"]').popover();

$('#geneButton').click(function(){
    $("input[type='file']").trigger('click');
})

// when model form submitted, show spinner (both for local run and remote trigger)
$("#model-form").on("submit", function () {
    $('#divLoading').show();
});

$('input:file').change(
    function(){
        if ($(this).val()) {
            $('#filename').text(this.value.replace(/C:\\fakepath\\/i, ''))
            $('input:submit').attr('disabled',false);
            $('#geneBtn').prop('disabled',true);
            console.log()
            var file = this.files[0];
            uploadFile(file);
    }
});

function uploadFile(file){
    var formData = new FormData();
    formData.append('formData', file);
    $.ajax({
        url: '/uploadgenes',  //Server script to process data
        type: 'POST',
        data: formData,
        contentType: false,
        processData: false
        //Ajax events
        //success: function(html){
        //    alert(html);
        //}
    });
    event.preventDefault();
}

function sampleGenes(){
    var your_array = ['CCNO','CENPF','LRRC56','ODAD3','DNAAF1','DNAAF6','DNAAF4','DNAH5','DNAH9','CFAP221','RSPH9',
    'FOXJ1','LRRC6','GAS2L2','DNAH1','GAS8','DNAI1','STK36','MCIDAS','RSPH4A','DNAAF3','DNAJB13','CCDC103','NME8',
    'ZMYND10','HYDIN','DNAAF5','CCDC40','ODAD2','DNAAF2','IFT122','INPP5E','CFAP298','DNAI2','SPAG1','SPEF2','ODAD4',
    'DNAL1','RSPH3','OFD1','CFAP300','CCDC65','DNAH11','RSPH1','DRC1','ODAD1'];
    var textarea = document.getElementById("enterGenes");
    textarea.value = your_array.join("\n");
}

function saveGenes(){
    var genes = $('#enterGenes').val();

    $.ajax({
        data: {
            genes: genes
        },
        type: 'POST',
        url: '/postgenes',
        success: function(data) {
            console.log(data);
            console.log(data.success);
            $('#geneModal').modal('hide');
            $('input:submit').attr('disabled',false);
            $("#geneButton").css("pointer-events", "none");
            $('#filename').empty();
        }
    });
    event.preventDefault();
}

function clearInput(){
    $.ajax({
        type: 'POST',
        url: '/clearinput',
        success: function(data) {
            console.log(data.success);
            $('input:submit').attr('disabled',true);
            $('#filename').val("");
            $('#geneBtn').prop('disabled', false);
            $("#geneButton").css("pointer-events", "auto");
        }
    });
    event.preventDefault();
}


function appendPrefix(){
    var jobname = $('#job').val();
    var prefix = $('#prefix').val();
    $.ajax({
        data: {
            jobname: jobname,
            prefix: prefix
        },
      type: 'POST',
      url: "/appendprefix",

      success: function(data){
        console.log(data.success);
        $('#job').val(data.jobname);
      }
    });
    event.preventDefault();
};

// run batch now works through same flask route as run local
// function runBatch(){
//     var jobname = $('#job').val();
//     var nettype = $('#network').val();
//     var features = $('#features').val();
//     var GSC = $('#negativeclass').val();

//     $.ajax({
//         data: {
//             jobname: jobname,
//             network: nettype,
//             feature: features,
//             negativeclass: GSC,
//         },
//       type: 'POST',
//       url: "/runbatch",
//       success: function(data){
//         console.log(data.success);
//       }
//     });
//     event.preventDefault();
// };