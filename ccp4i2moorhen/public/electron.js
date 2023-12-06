const path = require("path")
const express = require('express')
const process = require('process')
const http = require('http')
const qs = require('querystring')
const multer = require('multer')
const FormData = require("form-data")
const fileUpload = multer()
const { Readable } = require('stream')

const { app, BrowserWindow } = require("electron");
const isDev = require("electron-is-dev");
const { Console } = require("console");
const fetch = (...args) =>
	import('node-fetch').then(({default: fetch}) => fetch(...args));
// Handle creating/removing shortcuts on Windows when installing/uninstalling
if (require("electron-squirrel-startup")) {
    app.quit();
}

// Conditionally include the dev tools installer to load React Dev Tools
let installExtension, REACT_DEVELOPER_TOOLS;

if (isDev) {
    const devTools = require("electron-devtools-installer");
    installExtension = devTools.default;
    REACT_DEVELOPER_TOOLS = devTools.REACT_DEVELOPER_TOOLS;
}

function createWindow() {
    let CCP4i2Port = 43434
    let cootJob = null;
    // Create the browser window.
    const win = new BrowserWindow({
        width: 800,
        height: 600,
        icon: path.join(__dirname, "..", "src", "icons", "png", "128x128.png"),
        webPreferences: {
            nodeIntegration: true
        }
    });

    if (process.argv.length > 1) {
        CCP4i2Port = parseInt(process.argv[1])
    }

    if (process.argv.length > 2) {
        cootJob = process.argv[2]
    }

    let server;

    if (!isDev) {

        const MINPORT = 32778;
        const MAXPORT = 32800;

        const exp = express();

        //Install CORS headers
        exp.use(function (req, res, next) {
            res.header("Cross-Origin-Embedder-Policy", "require-corp");
            res.header("Cross-Origin-Opener-Policy", "same-origin");
            next();
        });

        exp.post("/database/uploadFileToJob", fileUpload.single('file'), async (oreq, ores) => {

            const formData = new FormData()
            Object.keys(oreq.body).forEach(key => {
                formData.append(key, oreq.body[key])
            })
            const fileStream = Readable.from(oreq.file.buffer);
            formData.append('file', fileStream)

            const newHeaders = {...oreq.headers, ...formData.getHeaders()}
            fetch(`http://127.0.0.1:${CCP4i2Port}/database/uploadFileToJob`, {
                method: "POST",
                body: formData,
                headers: newHeaders
            }).then(result => result.json())
                .then(jsonOfResult => {
                    ores.writeHead(result.statusCode, result.headers);
                    ores.write(JSON.stringify(jsonOfResult))
                    ores.end()
                })
                .catch(err => {
                    console.log(err)
                    console.log(`Error handling fetch with address ${"http://127.0.0.1:43434/database/uploadFileToJob"}`)
                })
        })

        //Forward database calls
        exp.get(/\/database|\/api/, (oreq, ores) => {
            console.log(`Handling database GET call`)
            const options = {
                host: '127.0.0.1', port: CCP4i2Port, path: oreq.url,
                method: oreq.method, headers: oreq.headers
            };

            const creq = http
                .request(options, pres => {
                    ores.writeHead(pres.statusCode, pres.headers);
                    pres.on('data', chunk => { ores.write(chunk); });
                    pres.on('close', () => { ores.end() });
                    pres.on('end', () => { ores.end() });
                })
                .on('error', e => {
                    console.log(e.message);
                    try {
                        ores.writeHead(500);
                        ores.write(e.message);
                    } catch (e) {
                        // ignore
                    }
                    ores.end();
                });
            creq.end();
        });

        //Is this something in the static tree ?
        exp.use(express.static(path.join(__dirname, "..", "build")));

        exp.get("/manifest.json/", (req, res) => {
            console.log(`mapping ${req.url} to manifest.json`)
            res.sendFile(path.join(__dirname, "..", "build", "manifest.json"))
        });

        exp.get(/\/*/, (req, res) => {
            console.log(`mapping ${req.url} to index`)
            res.sendFile(path.join(__dirname, "..", "build", "index.html"))
        });

        function serve(port) {
            server = exp.listen(port, () => {
                console.log('Listening on port:', server.address().port);
                win.loadURL(`http://localhost:${server.address().port}/?cootJob=${cootJob}&CCP4i2Port=${CCP4i2Port}`);
                //win.webContents.openDevTools()
            }).on('error', function (err) {
                if (port < MAXPORT) {
                    serve(port + 1);
                } else {
                    throw new Error("Run out of ports in Moorhen's range 32778-32800");
                }
            });
        }
        serve(MINPORT);
    } else {
        win.loadURL("http://localhost:9999");
    }


    // Open the DevTools.
    if (isDev) {
        win.webContents.openDevTools({ mode: "detach" });
    }
}

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.whenReady().then(createWindow);

// Quit when all windows are closed, except on macOS. There, it's common
// for applications and their menu bar to stay active until the user quits
// explicitly with Cmd + Q.
app.on("window-all-closed", () => {
    if (process.platform !== "darwin") {
        app.quit();
    }
});

app.on("activate", () => {
    // On macOS it's common to re-create a window in the app when the
    // dock icon is clicked and there are no other windows open.
    if (BrowserWindow.getAllWindows().length === 0) {
        createWindow();
    }
});

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.

