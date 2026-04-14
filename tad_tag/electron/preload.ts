import { contextBridge, ipcRenderer } from 'electron'

contextBridge.exposeInMainWorld('electronAPI', {
  // 1. 基础通信：监听主进程发来的信号 (包括更新进度和后端状态)
  on: (channel, callback) => {
    const validChannels = [
      'server-ready',
      'server-log',
      'main-process-message',
      'tad-finished',
      // --- 新增：自动更新相关的监听通道 ---
      'update-available',     // 发现新版本
      'update-not-available', // 已是最新版
      'update-progress',      // 下载进度百分比
      'update-downloaded',    // 下载完成
      'update-error'          // 更新出错
    ]

    if (validChannels.includes(channel)) {
      const subscription = (event, ...args) => callback(...args)
      ipcRenderer.on(channel, subscription)
      return () => ipcRenderer.removeListener(channel, subscription)
    }
  },

  // 2. 发送指令到主进程 (用于触发检查和开始下载)
  send: (channel, data) => {
    const validChannels = [
      'check-for-updates', // 触发检查
      'download-update'    // 触发下载
    ]
    if (validChannels.includes(channel)) {
      ipcRenderer.send(channel, data)
    }
  },

  // 3. 获取应用版本号 (从主进程的 package.json 中异步获取)
  getAppVersion: () => ipcRenderer.invoke('api:get-app-version'),

  // 4. 原有 TAD 业务接口
  runTadDetection: (params) => ipcRenderer.invoke('api:runTadDetection', params),
  checkServerStatus: () => ipcRenderer.invoke('api:check-server-status'),
  openFile: (options) => ipcRenderer.invoke('dialog:openFile', options),
  registerWithPS: (filePath) => ipcRenderer.invoke('api:registerWithPS', filePath),
  getTadStatus: (token) => ipcRenderer.invoke('api:getTadStatus', token),

  // 移除监听器
  removeListener: (channel) => {
    const validChannels = ['server-log', 'tad-finished', 'update-progress']
    if (validChannels.includes(channel)) {
      ipcRenderer.removeAllListeners(channel)
    }
  },
})