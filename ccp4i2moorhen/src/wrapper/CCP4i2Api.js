import $ from 'jquery'

const makeDbFilePromise = async (fileOfType, mimeType, arg) => {
    const fileDict = {
        filetypeid__filetypename: mimeType,
    }
    const dbFileIdNodes = $(fileOfType).find("dbFileId")
    if (dbFileIdNodes.length === 1 && $(dbFileIdNodes[0]).text().length > 0) {
        fileDict['fileid'] = $(dbFileIdNodes[0]).text()
    }
    const subTypeNodes = $(fileOfType).find("subType")
    if (subTypeNodes.length === 1 && $(subTypeNodes[0]).text().length > 0) {
        fileDict['filesubtype'] = parseInt($(subTypeNodes[0]).text())
    }
    const anotationNodes = $(fileOfType).find("annotation")
    if (anotationNodes.length === 1 && $(anotationNodes[0]).text().length > 0) {
        fileDict['annotation'] = $(anotationNodes[0]).text()
    }
    const fileUrl = `/getFileWithPredicate?${$.param({ fileid: fileDict.fileid })}`
    const annotation = fileDict.annotation ? fileDict.annotation : "A file"
    const subType = fileDict.filesubtype ? fileDict.filesubtype : 1
    return false//readFile(fileUrl, mimeType, annotation, subType, arg)
}

const makeParamsFilePromise = async (fileOfType, mimeType, arg) => {
    const fileDict = {}
    const projectNodes = $(fileOfType).find("project")
    if (projectNodes.length === 1 && $(projectNodes[0]).text().length > 0) {
        fileDict['projectId'] = $(projectNodes[0]).text()
    }
    const relPathNodes = $(fileOfType).find("relPath")
    if (relPathNodes.length === 1 && $(relPathNodes[0]).text().length > 0) {
        fileDict['relPath'] = $(relPathNodes[0]).text()
    }
    const baseNameNodes = $(fileOfType).find("baseName")
    if (baseNameNodes.length === 1 && $(baseNameNodes[0]).text().length > 0) {
        fileDict['baseName'] = $(baseNameNodes[0]).text()
    }
    let subType = 1
    const subTypeNodes = $(fileOfType).find("subType")
    if (subTypeNodes.length === 1 && $(subTypeNodes[0]).text().length > 0) {
        subType = $(subTypeNodes[0]).text()
    }
    let annotation = "A file"
    const annotationNodes = $(fileOfType).find("annotation")
    if (annotationNodes.length === 1 && $(annotationNodes[0]).text().length > 0) {
        annotation = $(annotationNodes[0]).text()
    }
    const fileUrl = `/getProjectFileData?${$.param(fileDict)}`
    return false//readFile(fileUrl, mimeType, annotation, subType, arg)
}

export class CCP4i2Api {
    constructor(urlRoot) {
        this.urlRoot = urlRoot
        this.postCommand = async (call, kws) => {
            const formData = new FormData()
            for (let key in kws) {
                formData.set(key, kws[key])
            }
            return fetch(`${this.urlRoot}/${call}`, {
                method: "POST",
                body: formData
            }).then(response => response.json())
        }
        this.uploadFileForJobObject = async (jobId, objectPath, fileBlob, fileName) => {
            const uploadResult = await this.postCommand("uploadFileForJobObject", {
                jobId, objectPath, fileName, file: fileBlob
            })
            return uploadResult
        }
        this.setJobParameterValue = async (jobId, objectPath, value) => {
            const rootTag = objectPath.split('.').at(-1)
            const argDict = {
                jobId: jobId,
                objectPath: objectPath,
                valueXMLText: `<${rootTag}>${value}</${rootTag}>`
            }
            console.log({ argDict })
            return this.postCommand("setJobParameterByXML", argDict)
        }
        this.readFileFromParamsXML = async () => {

        }
    }
}

