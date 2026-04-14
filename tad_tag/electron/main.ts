import { app, BrowserWindow, ipcMain, dialog } from 'electron'
import { createRequire } from 'node:module'
import { fileURLToPath } from 'node:url'
import path from 'node:path'
import fs from 'node:fs'
import { spawn, exec, ChildProcess } from 'node:child_process'
import { autoUpdater } from 'electron-updater'
import kill from 'tree-kill'
import { dialog } from 'electron'

const require = createRequire(import.meta.url)
const __dirname = path.dirname(fileURLToPath(import.meta.url))

// -------------------- 环境配置 --------------------
process.env.APP_ROOT = path.join(__dirname, '..')
export const VITE_DEV_SERVER_URL = process.env['VITE_DEV_SERVER_URL']
export const MAIN_DIST = path.join(process.env.APP_ROOT, 'dist-electron')
export const RENDERER_DIST = path.join(process.env.APP_ROOT, 'dist')
process.env.VITE_PUBLIC = VITE_DEV_SERVER_URL ? path.join(process.env.APP_ROOT, 'public') : RENDERER_DIST

let win: BrowserWindow | null = null
let serverProcess: ChildProcess | null = null
const pkg = require('../package.json')
const APP_VERSION = pkg.version

// 1. 定义全局变量记录状态
let isBackendReady = false;
// --- 工具函数：物理清洗路径（去除 Windows 不可见字符及统一斜杠） ---
const sanitizePath = (p: string) => {
  if (!p) return ''
  return p.replace(/[^\x20-\x7E]/g, '').replace(/\\/g, '/').trim()
}

// 3. 重要：添加 IPC 处理器，让前端能主动拿到 isBackendReady 的值
ipcMain.handle('api:check-server-status', () => {
  return isBackendReady;
});

// 在 createWindow 之后调用这个
function initUpdater(targetWin: BrowserWindow) {
  autoUpdater.autoDownload = false

  // 转发所有 updater 事件给渲染进程
  const sendStatus = (channel: string, data?: any) => {
    targetWin.webContents.send(channel, data)
  }

  autoUpdater.on('update-available', (info) => sendStatus('update-available', info))
  autoUpdater.on('update-not-available', () => sendStatus('update-not-available'))
  autoUpdater.on('download-progress', (p) => sendStatus('update-progress', p.percent))
  autoUpdater.on('update-downloaded', () => sendStatus('update-downloaded'))
  autoUpdater.on('error', (err) => sendStatus('update-error', err.message))
autoUpdater.on('update-downloaded', async (info) => {
    const dialogOpts = {
      type: 'info',
      buttons: ['立即重启并安装', '稍后提醒我'],
      title: '系统更新',
      message: `新版本 v${info.version} 已准备就绪`,
      detail: '更新包已下载完成。重启应用即可完成安装！',
      noLink: true
    }

    // 弹出原生对话框
    const { response } = await dialog.showMessageBox(targetWin, dialogOpts)

    if (response === 0) {
      // 用户点击了“立即重启”
      // 在安装前确保后端进程已关闭
      stopPackagedServer()
      autoUpdater.quitAndInstall()
    }
  })
}

// 响应 Preload 传来的指令
ipcMain.on('check-for-updates', () => autoUpdater.checkForUpdatesAndNotify())
ipcMain.on('download-update', () => autoUpdater.downloadUpdate())
// 将 'get-app-version' 修改为 'api:get-app-version' 以匹配 preload
ipcMain.handle('api:get-app-version', () => pkg.version)
// -------------------- Server 核心管理 --------------------

/**
 * 核心逻辑：轮询后端健康检查接口
 */
async function waitForServer(url: string, maxRetries = 20): Promise<boolean> {
  for (let i = 0; i < maxRetries; i++) {
    try {
      const response = await fetch(url)
      if (response.ok) return true
    } catch (e) {
      // 继续等待后端启动
    }
    await new Promise(resolve => setTimeout(resolve, 600))
  }
  return false
}

function startPackagedServer() {
  try {
    const serverConfig = {
      devPath: path.join(process.cwd(), 'backend'),
      prodPath: path.join(process.resourcesPath, 'backend'),
      filename: {
        win32: 'main.exe',
        darwin: 'main',
        linux: 'main'
      },
      args: []
    }

    const exeName = serverConfig.filename[process.platform as keyof typeof serverConfig.filename]
    const serverDir = app.isPackaged ? serverConfig.prodPath : serverConfig.devPath
    const serverExePath = path.join(serverDir, exeName)

    if (!fs.existsSync(serverExePath)) {
      console.error(`[Server] 可执行文件不存在: ${serverExePath}`)
      return false
    }

    if (process.platform !== 'win32') {
      try { fs.chmodSync(serverExePath, '755') } catch (err) {}
    }

    console.log(`[Server] 启动路径: ${serverExePath}`)

    // 优化：分离 stdio 避免主进程压力
    serverProcess = spawn(serverExePath, serverConfig.args, {
      shell: true, // 增强 Windows 稳定性
      stdio: ['ignore', 'pipe', 'pipe'],
      windowsHide: true,
      cwd: serverDir,
      env: { ...process.env,
    FORCE_COLOR: '1',
  PYTHONIOENCODING: 'utf-8',
      PYTHONUTF8: '1',
      LANG: 'en_US.UTF-8'
}
    })
    let logBuffer = '';

    serverProcess.stdout?.on('data', (data) => {
     const str = data.toString();
  // 实时发送给前端
  win?.webContents.send('server-log', str);
    })

    serverProcess.stderr?.on('data', (data) => {
      win?.webContents.send('server-log', ` ${data.toString().trim()}`)
    })

    serverProcess.on('close', (code) => {
      console.log(`[Server] 进程关闭，码: ${code}`)
      serverProcess = null
    })

    return true
  } catch (err: any) {
    console.error(`[Server] 启动异常:`, err)
    return false
  }
}

function stopPackagedServer() {
  console.log('[Server] 清理后端进程...')

  if (serverProcess && serverProcess.pid) {
    kill(serverProcess.pid, 'SIGKILL')
    serverProcess = null
  }

  // 异步执行强杀，防止关闭窗口时卡死
  const killCmd = process.platform === 'win32'
    ? 'taskkill /F /T /IM main.exe >nul 2>&1'
    : 'lsof -ti:5001 | xargs kill -9 >/dev/null 2>&1'

  exec(killCmd)
}

// -------------------- 窗口创建 --------------------

function createWindow() {
  win = new BrowserWindow({
    width: 1500,
    height: 1200,
    icon: path.join(process.env.VITE_PUBLIC!, 'TAD-TAG.ico'),
    webPreferences: {
      preload: path.join(__dirname, 'preload.mjs'),
      nodeIntegration: false,
      contextIsolation: true,
    },
  })
initUpdater(win)
  if (VITE_DEV_SERVER_URL) {
    win.loadURL(VITE_DEV_SERVER_URL)
    win.webContents.openDevTools()
  } else {
    win.loadFile(path.join(RENDERER_DIST, 'index.html'))
  }
}

// -------------------- 生命周期管理 --------------------

app.whenReady().then(async () => {
  createWindow()

  // 1. 启动后端
  startPackagedServer()

  const result = await waitForServer('http://127.0.0.1:5001/api/health');
  if (result) {
    isBackendReady = true;
      // 等页面加载完再发事件，避免事件丢失
    win?.webContents.once('did-finish-load', () => {
      win?.webContents.send('server-ready');
    });

    // 如果页面已经加载完了，did-finish-load 不会再触发，直接发
    if (!win?.webContents.isLoading()) {
      win?.webContents.send('server-ready');
    }
  }
})

app.on('window-all-closed', () => {
  stopPackagedServer()
  if (process.platform !== 'darwin') app.quit()
})








app.on('will-quit', () => {
  stopPackagedServer()
})

// -------------------- IPC 处理器 --------------------

// 通用的路径选择对话框
ipcMain.handle('dialog:openFile', async (event, options) => {
  const { canceled, filePaths } = await dialog.showOpenDialog(win!, options)
  if (canceled) return null
  return filePaths.map(p => sanitizePath(p))
})

// 使用 PowerShell 转发请求（应对复杂路径或编码问题）
ipcMain.handle('api:registerWithPS', async (event, filePath) => {
  return new Promise((resolve) => {
    const cleanPath = sanitizePath(filePath)

    // 构造 PowerShell 原生调用
    const cmd = `powershell -Command "Invoke-RestMethod -Uri 'http://127.0.0.1:5001/api/register' -Method Post -ContentType 'application/json; charset=utf-8' -Body '{\\\"path\\\": \\\"${cleanPath}\\\", \\\"type\\\": \\\"hic\\\"}'"`

    console.log(`[Executing PS]: ${cmd}`)

    exec(cmd, (error, stdout, stderr) => {
      if (error) {
        resolve({ success: false, message: stderr || error.message })
        return
      }
      try {
        const data = JSON.parse(stdout)
        resolve({ success: true, data })
      } catch (e) {
        resolve({ success: false, message: 'JSON 解析失败', raw: stdout })
      }
    })
  })
})

// 运行 TAD 检测
//ipcMain.handle('api:runTadDetection', async (event, tadParams) => {
// try {
//   if (tadParams.path) tadParams.path = sanitizePath(tadParams.path)
//
//   const response = await fetch('http://127.0.0.1:5001/api/tad/run', {
//     method: 'POST',
//     headers: { 'Content-Type': 'application/json' },
//     body: JSON.stringify(tadParams)
//   })
//   const data = await response.json()
//   return { success: response.ok, data }
// } catch (err: any) {
//   return { success: false, message: `Fetch Error: ${err.message}` }
// }
//})

ipcMain.handle('api:runTadDetection', async (event, tadParams) => {
  try {
    if (tadParams.path) tadParams.path = sanitizePath(tadParams.path);

    // 1. 发送启动指令
    const response = await fetch('http://127.0.0.1:5001/api/tad/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(tadParams)
    });

    const startData = await response.json();

    if (!response.ok) throw new Error(startData.error || 'Failed to start task');

    // 2. 既然任务已经开始，主进程开启一个后台监控逻辑
    // 这样渲染进程（前端）就不需要一直维持那个 HTTP 连接
    monitorTadStatus(event, tadParams.token);

    // 3. 立即给前端返回“任务已启动”的消息
    return { success: true, message: 'Task started in background', data: startData };

  } catch (err) {
    return { success: false, message: `Launch Error: ${err.message}` };
  }
});

// 后台监控函数
async function monitorTadStatus(event, token) {
  const check = async () => {
    try {
      const res = await fetch(`http://127.0.0.1:5001/api/tad/status/${token}`);
      const statusData = await res.json();

      if (statusData.status === 'done') {
        // 任务成功，通过之前监听的 server-log 或新频道发送通知
        event.sender.send('tad-finished', { success: true, token });
      } else if (statusData.status === 'failed') {
        event.sender.send('tad-finished', { success: false, error: statusData.error_msg });
      } else {
        // 还在运行，继续轮询（例如每 10 秒检查一次）
        setTimeout(check, 10000);
      }
    } catch (e) {
      console.error('Polling error:', e);
      setTimeout(check, 20000); // 出错则放慢频率
    }
  };
  check();
}

// 获取任务状态
ipcMain.handle('api:getTadStatus', async (event, token) => {
  try {
    const response = await fetch(`http://127.0.0.1:5001/api/tad/status/${token}`)
    const data = await response.json()
    return { success: response.ok, data }
  } catch (err: any) {
    return { success: false, message: err.message }
  }
})