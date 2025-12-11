var queryStack = new Array()
var queryPointer = 0;
var forwardDiv = null
var backwardDiv = null
var lineColours = ['red','green','blue','cyan','magenta','yellow','black']

function onLoad(){
    //Throw away all existing elements
    bodyDiv =document.getElementsByTagName('body')[0]
    while (bodyDiv.firstChild){
        bodyDiv.removeChild(bodyDiv.firstChild)
    }
    
    bodyDiv.innerHTML="\
    <div id='bannerDiv' style='height:100px;background-color:red;'> \
    <form id='file-form' action='handler.php' method='POST'>\
    <input type='file' id='file-select' name='zippedProjects[]'/>\
    <button type='submit' id='upload-button'>Upload</button><progress style='width:250px;'/>\
    </form><div id='feedbackDiv' style='height:50px; width:200px;'></div></div>\
    <div id='navigationDiv' style ='height:50px'> </div>\
    <div id='treeBar' style='width:400px;float:left;overflow:auto;'> </div>\
    <div id='contentDiv' style='float:left;background-cyan;overflow:scroll;width:100%;'> </div>\
    "
    
    // Get the selected files from the input.
    var form = document.getElementById('file-form');
    var fileSelect = document.getElementById('file-select');
    var uploadButton = document.getElementById('upload-button');
    form.onsubmit = function(event) {
        event.preventDefault();
        
        // Update button text.
        uploadButton.innerHTML = 'Uploading...';
        
        // The rest of the code will go here...
        var files = fileSelect.files;
        // Create a new FormData object.
        var formData = new FormData();
        
        // Loop through each of the selected files.
        for (var i = 0; i < 1; i++) {
            var file = files[i];
            
            // Add the file to the request.
            formData.append('file', file, file.name);
        }
        // Set up the request.
        var xhr = new XMLHttpRequest();
        xhr.upload.addEventListener('progress',progressHandlingFunction, false);
        // Open the connection.
        xhr.open('POST', 'handler.php', true);
        // Set up a handler for when the request finishes.
        xhr.onload = function () {
            if (xhr.status === 200) {
                // File(s) uploaded.
                uploadButton.innerHTML = 'Upload';
                ajaxRetrieve("?listProjects", projectListRetrieved)
            } else {
                alert('An error occurred!');
            }
        };
        console.log('xhr',xhr)
        console.log('formData',formData);
        // Send the Data.
        xhr.send(formData);
    }
    
    
    var contentDiv = document.getElementById('contentDiv')
    ajaxRetrieve("?listProjects", projectListRetrieved)
    window.onresize = function(){
        onWindowResize()
    }
    onWindowResize()
}

function progressHandlingFunction(e){
    if(e.lengthComputable){
        $('progress').attr({value:e.loaded,max:e.total});
    }
}
function clearContentDiv() {
    var contentDiv = document.getElementById('contentDiv')
    //Remove existing content from contentDiv
    while (contentDiv.firstChild) {
        contentDiv.removeChild(contentDiv.firstChild);
    }
}

function updateNavigationDiv(){
    /*Remove existing content from contentDiv
     while (navigationDiv.firstChild) {
     navigationDiv.removeChild(navigationDiv.firstChild);
     }*/
    //forwardDiv.removeChild(forwardDiv.firstChild)
    //backwardDiv.removeChild(backwardDiv.firstChild)
    var navigationDiv = document.getElementById('navigationDiv')
    if (forwardDiv == null){
        forwardDiv = document.createElement('div');
        forwardDiv.style.cssText="width:80px; height:40px; float:left;";
        navigationDiv.appendChild(forwardDiv);
    }
    else {
        while (forwardDiv.firstChild){ forwardDiv.removeChild(forwardDiv.firstChild) }
    }
    if (backwardDiv == null){
        backwardDiv = document.createElement('div');
        backwardDiv.style.cssText="width:80px; height:40px; float:left;";
        navigationDiv.appendChild(backwardDiv);
    }
    else {
        while (backwardDiv.firstChild){ backwardDiv.removeChild(backwardDiv.firstChild) }
    }
    if (queryStack.length > queryPointer) {
        forwardButton = document.createElement('BUTTON')
        forwardDiv.appendChild(forwardButton)
        forwardButtonText = document.createTextNode('Forward')
        forwardButton.appendChild(forwardButtonText)
        forwardButton.onclick = forward
    }
    if (queryPointer > 1) {
        backButton = document.createElement('BUTTON')
        backwardDiv.appendChild(backButton)
        backButtonText = document.createTextNode('Back')
        backButton.appendChild(backButtonText)
        backButton.onclick = back
    }
}

function projectListRetrieved(responseText){
    //alert(responseText)
    var projectTupleArray = JSON.parse(responseText);
    projectDicts = []
    for (var iProject = 0; iProject < projectTupleArray.length; iProject++){
        projectTuple = projectTupleArray[iProject];
        projectDicts.push({ProjectID:projectTuple[0], ProjectName:projectTuple[1], ProjectDirectory:projectTuple[2],ParentID:projectTuple[3]});
    }
    
    //Empty the content Div
    clearContentDiv()
    //Add navigation buttons as appropriate
    updateNavigationDiv()
    
    var idToNameAssociativeArray = {}
    for (var projectDictNumber=0; projectDictNumber<projectDicts.length; projectDictNumber++){
        var projectDict = projectDicts[projectDictNumber]
        idToNameAssociativeArray[projectDict['ProjectID']] = projectDict['ProjectName']
        
        projectDict['id'] = projectDict['ProjectID']
        projectDict['text'] = projectDict['ProjectName']
        projectDict['state'] = {'opened':true}
        if (projectDict['ParentID'] != null){
            projectDict['parent'] = projectDict['ParentID']
        }
        else {
            projectDict['parent'] = '#'
        }
        
    }
    var contentDiv = document.getElementById('contentDiv')
    var table=document.createElement('table');
    contentDiv.appendChild(table)
    table.style.width='100%';
    table.setAttribute('border','1');
    var tbody=document.createElement('tbody');
    tbody.class='fancy'
    table.appendChild(tbody)
    
    var tr=document.createElement('tr');
    tbody.appendChild(tr)
    var headersAndTags = [['Project name','ProjectName'],['ProjectID','ProjectID'],['Parent project name','ParentID'],['Get jobs','Get jobs'] ,['Get files','Get files']]
    for (var iHeader=0; iHeader<headersAndTags.length; iHeader++){
        var th=document.createElement('th');
        tr.appendChild(th)
        th.innerHTML = headersAndTags[iHeader][0]
    }
    
    for (var projectDictNumber=0; projectDictNumber<projectDicts.length; projectDictNumber++){
        projectDict = projectDicts[projectDictNumber]
        
        var tr=document.createElement('tr');
        tbody.appendChild(tr)
        
        
        for (var iHeader=0; iHeader<headersAndTags.length; iHeader++){
            header = headersAndTags[iHeader][0];
            tag = headersAndTags[iHeader][1]
            var td=document.createElement('td');
            tr.appendChild(td)
            if (tag == 'ParentID'){
                text = document.createTextNode(idToNameAssociativeArray[projectDict[tag]])
                td.appendChild(text)
            }
            else if (tag == 'Get jobs'){
                button = document.createElement('BUTTON')
                button.setAttribute('id', projectDict['ProjectID'])
                button.onclick = getProjectJobsClicked
                td.appendChild(button)
                buttonText = document.createTextNode('Fetch files')
                button.appendChild(buttonText)
            }
            else if (tag == 'Get files'){
                button = document.createElement('BUTTON')
                button.setAttribute('id', projectDict['ProjectID'])
                button.onclick = getProjectFilesClicked
                td.appendChild(button)
                buttonText = document.createTextNode('Fetch files')
                button.appendChild(buttonText)
            }
            else {
                text = document.createTextNode(projectDict[tag])
                td.appendChild(text)
            }
        }
        
    }
    //Try to destroy tree in the treeBar. Wrap attempt try...catch because there may not yet be one !
    try{
        $('#treeBar').jstree('destroy')
    }
    catch(err){
    }
    //alert(JSON.stringify(projectDicts))
    $('#treeBar')// listen for event
    .on('changed.jstree', function (e, data) {
        var i, j, r = [];
        for(i = 0, j = data.selected.length; i < j; i++) {
        node =data.instance.get_node(data.selected[i])
        r.push(JSON.stringify(node.original.text));
        fakeEvent = new Object()
        fakeEvent.target = new Object()
        fakeEvent.target.id = node.original.ProjectID
        getProjectJobsClicked(fakeEvent)
        }
        $('#feedbackDiv').html('Selected: ' + r.join(', '));
        }).jstree({ 'core' : {
                  'data' : projectDicts,
                  'check_callback' : true
                  },
                  'plugins':['crrm','dnd']
                  }).bind("move_node.jstree", function (e, data) {
                          ajaxQuery="?updateProject?projectId="+data.node.id+"?key=parentprojectid?value="+data.parent;
                          ajaxRetrieve(ajaxQuery);
                          });
}

function clearContentDiv() {
    var contentDiv = document.getElementById('contentDiv')
    //Remove existing content from contentDiv
    while (contentDiv.firstChild) {
        contentDiv.removeChild(contentDiv.firstChild);
    }
}

function getProjectJobsClicked(event) {
    projectId = event.target.id
    ajaxRetrieve("?getProjectJobListInfo?topLevelOnly=False?projectId="+projectId, projectJobsRetrieved)
}

function projectJobsRetrieved(responseText){
    jobInfoList = JSON.parse(responseText)
    
    for (var jobNumber=0; jobNumber<jobInfoList.length; jobNumber++){
        var jobDict = jobInfoList[jobNumber]
        jobDict['id'] = jobDict['jobid']
        if (jobDict['jobtitle'] != null){
            jobDict['text'] = jobDict['jobnumber']+jobDict['jobtitle']
        }
        else {
            jobDict['text'] = jobDict['jobnumber']+jobDict['taskname']
        }
        jobDict['state'] = {'opened':true}
        if (jobDict['status']=='Finished'){
            jobDict['state'].disabled=false
        }
        else {
            jobDict['state'].disabled=true
        }
        if (jobDict['parentjobid'] != null){
            jobDict['parent'] = jobDict['parentjobid']
        }
        else {
            jobDict['parent'] = '#'
            jobDict['state'].opened=false
        }
    }
    
    //Empty the content Div
    clearContentDiv()
    //Add navigation buttons as appropriate
    updateNavigationDiv()
    //Create divs for table and content
    var tableDiv = document.createElement('div')
    tableDiv.setAttribute('id', 'tableDiv');
    tableDiv.style.cssText = 'overflow:scroll;border-width:1px;border-style:solid;border-color:black;height:100%;'
    contentDiv.appendChild(tableDiv)
    
    var table=document.createElement('table');
    tableDiv.appendChild(table)
    table.style.width='100%';
    table.setAttribute('border','1');
    var tbody=document.createElement('tbody');
    tbody.class='fancy'
    table.appendChild(tbody)
    var tr=document.createElement('tr');
    tbody.appendChild(tr)
    var headersAndKeys = [
                          ['Job Number','jobnumber'],
                          ['Title','jobtitle'],
                          ['Task name','taskname'],
                          ['Status','status'],
                          ['Created at','creationtime'],
                          ['Report','report']
                          ]
    for (var iHeader=0; iHeader<headersAndKeys.length; iHeader++){
        var th=document.createElement('th');
        tr.appendChild(th)
        text = document.createTextNode(headersAndKeys[iHeader][0])
        th.appendChild(text)
    }
    for (var iJobInfo=0; iJobInfo<jobInfoList.length; iJobInfo++){
        jobInfo = jobInfoList[iJobInfo]
        var tr=document.createElement('tr');
        tbody.appendChild(tr)
        for (var iHeader=0; iHeader<headersAndKeys.length-1; iHeader++){
            td=document.createElement('td');
            tr.appendChild(td)
            text = document.createTextNode(jobInfo[headersAndKeys[iHeader][1]])
            td.appendChild(text)
        }
        td = document.createElement('td')
        tr.appendChild(td)
        if (jobInfo['status'] === 'Finished'){
            button = document.createElement('BUTTON')
            button.setAttribute('id', jobInfo['jobid'])
            button.onclick = getReportClicked
            td.appendChild(button)
            buttonText = document.createTextNode('Fetch report')
            button.appendChild(buttonText)
        }
    }
    //[{"finishtime": null, "status": "Pending", "childjobs": [], "preceedingjobid": null, "taskversion": null, "parentjobid": null, "jobtitle": "Coot execute script", "creationtime": "09:41 13/Dec/13", "jobid": "aee626ba63da11e3bc1eb8e85634974a", "projectid": "76918807fdad11e2a09358b0357c3c6c", "useragent": "CCP4i2", "taskname": "coot_script_lines", "evaluation": "Unknown", "jobnumber": "49"}]
    //alert(JSON.stringify(projectDicts))
    //treeBar = document.getElementById('treeBar')
    $('#treeBar').jstree('destroy')
    $('#treeBar').on('changed.jstree', function (e, data) {
                     var i, j, r = [];
                     for(i = 0, j = data.selected.length; i < j; i++) {
                     var node =data.instance.get_node(data.selected[i])
                     r.push(JSON.stringify(node.original.text));
                     var fakeEvent = {target:{id:node.original.jobid}};
                     getReportClicked(fakeEvent)
                     }
                     $('#feedbackDiv').html('Selected: ' + r.join(', '));
                     }).jstree({ 'core' : {
                               'data' : jobInfoList
                               } });
    
}

function iFrameGetReportClicked(event){
    jobId = event.target.id;
    var contentDiv = document.getElementById('contentDiv');
    while (contentDiv.firstChild){
        contentDiv.removeChild(contentDiv.firstChild);
    }
    iframe = document.createElement("object");
    iframe.style.height="100%";
    iframe.style.width="100%";
    iframe.setAttribute("src", "/database"+"?File?jobId="+jobId+"?filePath=report.html");
    contentDiv.appendChild(iframe);
}

function getReportClicked(event) {
    jobId = event.target.id
    ajaxRetrieve("?File?jobId="+jobId+"?filePath=report.html", jobReportRetrieved)
}

function jobReportRetrieved(responseXML){
    //alert(new XMLSerializer().serializeToString(responseXML))
    //Empty the content Div
    clearContentDiv()
    //Andmove across cloned elements from the body element of the responseXML
    var contentDiv = document.getElementById('contentDiv')
    
    var $xml = $( responseXML );
    var bodyNode = $xml.find( "body" ).get(0);
    
    for (var iChild=0; iChild< bodyNode.childNodes.length; iChild++){
        child=bodyNode.childNodes[iChild];
        contentDiv.innerHTML += (new XMLSerializer().serializeToString(child))
    }
    
    console.log('Patching imgs');
    //Patch up img elements to fetch image from the Webapp part of the server
    $("img").each(function(){
                  var img = this;
                  var oldSrc = img.getAttribute('src');
                  //img.setAttribute('src','database?File?jobId='+jobId+'?filePath='+oldSrc);
                  });
    
    console.log('Patching drawn divs');
    //Patch up DrawnDiv elements that provide data as URLs so that those URLs reference the database
    $(":div[data-renderer]").each(function(){
                  var drawnDiv = this;
                  if ($(drawnDiv).data('is-urls') === 'True'){
                      oldUrl = $(drawnDiv).data('data');
                      newUrl = "/database"+"?File?jobId="+jobId+"?filePath="+oldUrl;
                      $(drawnDiv).get(0).setAttribute('data-data', newUrl);
                  }
                  console.log('About to draw drawn div...'+newUrl);
                                  try {drawDrawnDiv(drawnDiv);}
                                  catch(err) {console.log('Failed');}
                                  
                  console.log('Done');
                  });
    console.log('Done with drawn divs');
    
    //Patch up file download objects
    $("object.qt_object").each(function(){
                               var params = {};
                               for (iChildNode = 0; iChildNode < this.childNodes.length; iChildNode++){
                               var childNode = this.childNodes[iChildNode];
                               params[childNode.getAttribute('name')] = childNode.getAttribute('value');
                               }
                               //alert (JSON.stringify(params))
                               newDiv = document.createElement("div")
                               newDiv.setAttribute('class','qt_object');
                               this.parentNode.replaceChild(newDiv, this);
                               newDiv.innerHTML = '<img src="/database?File?icon='+params['mimeTypeName']+'" height=20 width=20 style="margin:5px;display:inline-block;border:1px solid grey;"></img><div style="width:500px;height:25px;float:right;font-size:125%;border:0px solid red;margin:2px;" id="Download_'+params['dbFileId']+'" class="downloadable">'+params['annotation']+'</div>';
                               //newDiv.innerHTML = '<img src="/database?File?icon='+params['mimeTypeName']+'" height=20 width=20 style="margin:5px;display:inline-block;border:1px solid grey;"></img><span style="margin:5px;font-size:125%;" id="Download_'+params['dbFileId']+'">'+params['annotation']+'</span>';
                               var fileSpec = "/database?File?fileId="+params['dbFileId'];
                               $('#Download_'+params['dbFileId']).click(function(){
                                                                        window.downloadFile(fileSpec)
                                                                        });
                               newDiv.style.height = this.getAttribute("height")+'px';
                               newDiv.style.width = this.getAttribute("width")+'px';
                               newDiv.style.border = "1px solid green";
                               });
    
    //Patch up file download objects
    $("object.qt_launch").each(function(){
                               if (this.getAttribute('type') === 'x-ccp4-widget/CDownloadButton'){
                               var params = {};
                               for (iChildNode = 0; iChildNode < this.childNodes.length; iChildNode++){
                               var childNode = this.childNodes[iChildNode];
                               params[childNode.getAttribute('name')] = childNode.getAttribute('value');
                               }
                               if ('jobId' in params && 'dataName' in params){
                               var fileSpec = "/database?File?jobId="+params['jobId']+"?filePath=tables_as_csv_files/"+params['dataName']+".csv";
                               var newButton =document.createElement('Button');
                               newButton.innerHTML = 'Download';
                               $(newButton).click(function(){
                                                  window.downloadFile(fileSpec)
                                                  });
                               this.parentNode.replaceChild(newButton, this);
                               }
                               }
                               });
    
    //Add navigation buttons as appropriate
    updateNavigationDiv()
    // For script elements, hunt to see if one is an entry into require js
    var scriptElements = responseXML.getElementsByTagName('script')
    for (var iScript=0; iScript < scriptElements.length; iScript++){
        var scriptElement = scriptElements[iScript];
        if (scriptElement.hasAttribute('data-main')){
            var requireRootBits = scriptElement.getAttribute('data-main').split('/');
            var lastBit = requireRootBits[requireRootBits.length-1]
            require([lastBit],function(){});
        }
    }
}

window.downloadFile = function(sUrl) {
    
    //If in Chrome or Safari - download via virtual link click
    if (window.downloadFile.isChrome || window.downloadFile.isFirefox) {
        //Creating new link node.
        var link = document.createElement('a');
        link.href = sUrl;
        if (link.download !== undefined){
            //Set HTML5 download attribute. This will prevent file from opening if supported.
            var fileName = sUrl.substring(sUrl.lastIndexOf('/') + 1, sUrl.length);
            link.download = fileName;
        }
        
        //Dispatching click event.
        if (document.createEvent) {
            var e = document.createEvent('MouseEvents');
            e.initEvent('click' ,true ,true);
            link.dispatchEvent(e);
            return true;
        }
    }
    
    // Force file download (whether supported by server).
    var query = '?Download';
    
    window.location.assign(sUrl + query);
}
window.downloadFile.isChrome = navigator.userAgent.toLowerCase().indexOf('chrome') > -1;
window.downloadFile.isSafari = navigator.userAgent.toLowerCase().indexOf('safari') > -1;
window.downloadFile.isFirefox = navigator.userAgent.toLowerCase().indexOf('firefox') > -1;

function getProjectFilesClicked(event) {
    projectId = event.target.id
    ajaxRetrieve("?getProjectFiles?projectId="+projectId, projectFilesRetrieved)
}

function projectFilesRetrieved(responseText){
    fileInfoList = JSON.parse(responseText)
    //0:Files.JobId, 1:Files.FileID, 2:ImportFiles.ImportId, 3:Files.FileTypeId, 4:Files.Filename, 5:Files.Annotation
    //Empty the content Div
    clearContentDiv()
    //Add navigation buttons as appropriate
    updateNavigationDiv()
    //Create divs for table and content
    var tableDiv = document.createElement('div')
    tableDiv.setAttribute('id', 'tableDiv');
    tableDiv.style.cssText = 'overflow:scroll;border-width:1px;border-style:solid;border-color:black;height:75%;'
    contentDiv.appendChild(tableDiv)
    var detailDiv = document.createElement('div')
    detailDiv.setAttribute('id', 'detailDiv')
    detailDiv.style.cssText = 'overflow:scroll;border-width:1px;border-style:solid;border-color:black;height:25%;'
    contentDiv.appendChild(detailDiv)
    
    var table=document.createElement('table');
    tableDiv.appendChild(table)
    table.style.width='100%';
    table.setAttribute('border','1');
    var tbody=document.createElement('tbody');
    tbody.class='fancy'
    table.appendChild(tbody)
    var tr=document.createElement('tr');
    tbody.appendChild(tr)
    var headersAndKeys = [
                          ['Job ID',0],
                          ['Filename',4],
                          ['Annotation',5],
                          ['File class',3]
                          ]
    for (var iHeader=0; iHeader < headersAndKeys.length; iHeader++){
        var th=document.createElement('th');
        tr.appendChild(th)
        text = document.createTextNode(headersAndKeys[iHeader][0])
        th.appendChild(text)
    }
    for (var iFileInfo=0; iFileInfo<fileInfoList.length; iFileInfo++){
        fileInfo = fileInfoList[iFileInfo]
        var tr=document.createElement('tr');
        tbody.appendChild(tr)
        for (var iHeader=0; iHeader<headersAndKeys.length; iHeader++){
            td=document.createElement('td');
            tr.appendChild(td)
            text = document.createTextNode(fileInfo[headersAndKeys[iHeader][1]])
            td.appendChild(text)
        }
    }
    //{"FileTypeId": 2, "JobID": "5b358a75000711e3be8f58b035fea8fd", "Filename": "XYZOUT.pdb", "ImportId": null, "FileID": "e977e2f0000711e3b02558b035fea8fd", "Annotation": null, "FileTypeClass": "PdbDataFile"}
}

function ajaxRetrieve(queryString, responseProcessor){
    //Increment query Pointer
    if (queryPointer+1>queryStack.length) {
        queryStack.push({'Query':queryString,'ResponseProcessor':responseProcessor})
        queryPointer += 1
    }
    else {
        oldQueryAtPosition = queryStack[queryPointer]
        if (queryString != oldQueryAtPosition['Query'] || responseProcessor != oldQueryAtPosition['ResponseProcessor']){
            queryStack = queryStack.slice(0,queryPointer)
            queryStack.push({'Query':queryString,'ResponseProcessor':responseProcessor})
        }
        queryPointer += 1
    }
    var mygetrequest = new ajaxRequest()
    mygetrequest.onreadystatechange=function(){
        if (mygetrequest.readyState==4){
            if (mygetrequest.status==200 || window.location.href.indexOf("http")==-1){
                if (mygetrequest.responseXML != null) {
                    responseProcessor(mygetrequest.responseXML)
                }
                else {
                    responseProcessor(mygetrequest.responseText)
                }
                
            }
        }
    }
    mygetrequest.open("GET", queryString, true)
    mygetrequest.send()
}

function forward()
{
    query = queryStack[queryPointer]
    ajaxRetrieve(query['Query'], query['ResponseProcessor'])
}

function back()
{
    query = queryStack[queryPointer-2]
    queryPointer -= 2
    ajaxRetrieve(query['Query'], query['ResponseProcessor'])
}

function ajaxRequest(){
    var activexmodes=["Msxml2.XMLHTTP", "Microsoft.XMLHTTP"] //activeX versions to check for in IE
    if (window.ActiveXObject){ //Test for support for ActiveXObject in IE first (as XMLHttpRequest in IE7 is broken)
        for (var i=0; i<activexmodes.length; i++){
            try{
                return new ActiveXObject(activexmodes[i])
            }
            catch(e){
                //suppress error
            }
        }
    }
    else if (window.XMLHttpRequest){ // if Mozilla, Safari etc
        return new XMLHttpRequest();
    }
    else {
        return false
    }
}


function onWindowResize(){
    document.getElementById("bannerDiv").style.width = (window.innerWidth-50)+"px"
    document.getElementById("navigationDiv").style.width = (window.innerWidth-50)+"px"
    document.getElementById("treeBar").style.height = (window.innerHeight - 200)+"px"
    document.getElementById("contentDiv").style.height = (window.innerHeight - 200)+"px"
    document.getElementById("contentDiv").style.width = (window.innerWidth - 440)+"px"
}

