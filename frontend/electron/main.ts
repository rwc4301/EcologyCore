import { app, BrowserWindow } from 'electron'
import * as path from 'path'
import * as isDev from 'electron-is-dev'

function createWindow() {
  const win = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false
    }
  })

  if (isDev) {
    win.loadURL('http://localhost:3000')
  } else {
    win.loadFile(path.join(__dirname, '../out/index.html'))
  }
}

app.whenReady().then(createWindow)

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit()
  }
})
