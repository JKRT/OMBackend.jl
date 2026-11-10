# The state of a stalled process, from outside (branch diag/windows-ext-hang only).
#   powershell -File snapshot.ps1 <prefix> <pid>
# Writes <prefix>.snapshot.txt (processes, CPU samples, thread wait reasons, loaded gettext/Glib libraries, window
# titles, popup events of the system log) and <prefix>.screen.png. Every section stands alone: one that fails
# says so and the next one runs.
param([string]$Prefix, [int]$TargetPid)

$ErrorActionPreference = 'Stop'
$out = "$Prefix.snapshot.txt"
function Section([string]$title, [scriptblock]$body) {
  "==== $title ====" | Out-File -FilePath $out -Append -Encoding utf8
  try {
    $text = (& $body | Out-String -Width 400)
    $text | Out-File -FilePath $out -Append -Encoding utf8
  } catch {
    "FAILED: $($_.Exception.Message)" | Out-File -FilePath $out -Append -Encoding utf8
  }
}

"snapshot of pid $TargetPid at $(Get-Date -Format o)" | Out-File -FilePath $out -Encoding utf8

Section "julia, WerFault and GR processes" {
  Get-CimInstance Win32_Process |
    Where-Object { $_.Name -match 'julia|WerFault|gksqt|csrss' } |
    Select-Object ProcessId, ParentProcessId, Name, ThreadCount, KernelModeTime, UserModeTime, CreationDate, CommandLine |
    Format-List
}

Section "CPU time of the process, two samples 5 s apart (no change: it waits)" {
  $a = (Get-Process -Id $TargetPid).TotalProcessorTime
  Start-Sleep -Seconds 5
  $b = (Get-Process -Id $TargetPid).TotalProcessorTime
  "first $a, second $b, used $(($b - $a).TotalMilliseconds) ms"
}

Section "threads of the process (a hard-error popup: a thread waiting for an LPC reply)" {
  (Get-Process -Id $TargetPid).Threads | ForEach-Object {
    $reason = ''
    try { if ("$($_.ThreadState)" -eq 'Wait') { $reason = "$($_.WaitReason)" } } catch { $reason = 'not readable' }
    "thread $($_.Id) state $($_.ThreadState) wait reason $reason"
  }
}

Section "gettext, iconv, Glib, GR and Modelica runtime libraries in the process" {
  (Get-Process -Id $TargetPid).Modules |
    Where-Object { $_.ModuleName -match 'intl|iconv|gio-|glib-|gobject-|gmodule-|Modelica|OpenModelica|libGR|gks|omcgc' } |
    Select-Object ModuleName, FileName |
    Format-Table -AutoSize
}

Section "the last libraries the process loaded (the end of its module list)" {
  (Get-Process -Id $TargetPid).Modules | Select-Object -Last 25 -Property ModuleName, FileName | Format-Table -AutoSize
}

Section "tasklist /v of julia.exe (window titles)" {
  tasklist /v /fi "imagename eq julia.exe"
}

Section "visible windows with a title (a dialog has the class #32770)" {
  Add-Type -TypeDefinition @'
using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public static class DiagWindows {
  public delegate bool EnumProc(IntPtr hWnd, IntPtr lParam);
  [DllImport("user32.dll")] static extern bool EnumWindows(EnumProc proc, IntPtr lParam);
  [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr hWnd);
  [DllImport("user32.dll", CharSet = CharSet.Unicode)] static extern int GetWindowText(IntPtr hWnd, StringBuilder text, int count);
  [DllImport("user32.dll", CharSet = CharSet.Unicode)] static extern int GetClassName(IntPtr hWnd, StringBuilder text, int count);
  [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr hWnd, out uint pid);
  public static List<string> List() {
    List<string> found = new List<string>();
    EnumWindows(delegate(IntPtr hWnd, IntPtr lParam) {
      StringBuilder title = new StringBuilder(512);
      StringBuilder cls = new StringBuilder(256);
      GetWindowText(hWnd, title, title.Capacity);
      GetClassName(hWnd, cls, cls.Capacity);
      uint pid;
      GetWindowThreadProcessId(hWnd, out pid);
      if (title.Length > 0 || cls.ToString() == "#32770") {
        found.Add("pid " + pid + " visible " + IsWindowVisible(hWnd) + " class " + cls + " title " + title);
      }
      return true;
    }, IntPtr.Zero);
    return found;
  }
}
'@
  [DiagWindows]::List()
}

Section "popups in the system log (Application Popup, event 26: the text of a hard-error dialog)" {
  Get-WinEvent -FilterHashtable @{ LogName = 'System'; Id = 26 } -MaxEvents 10 -ErrorAction Stop |
    Select-Object TimeCreated, ProviderName, Message |
    Format-List
}

Section "application errors in the application log (the last hour)" {
  Get-WinEvent -FilterHashtable @{ LogName = 'Application'; Level = 2; StartTime = (Get-Date).AddHours(-1) } -MaxEvents 15 -ErrorAction Stop |
    Select-Object TimeCreated, ProviderName, Id, Message |
    Format-List
}

Section "screenshot" {
  Add-Type -AssemblyName System.Windows.Forms
  Add-Type -AssemblyName System.Drawing
  $screen = [System.Windows.Forms.SystemInformation]::VirtualScreen
  $bitmap = New-Object System.Drawing.Bitmap $screen.Width, $screen.Height
  $graphics = [System.Drawing.Graphics]::FromImage($bitmap)
  $graphics.CopyFromScreen($screen.Left, $screen.Top, 0, 0, $bitmap.Size)
  $bitmap.Save("$Prefix.screen.png", [System.Drawing.Imaging.ImageFormat]::Png)
  $graphics.Dispose()
  $bitmap.Dispose()
  "saved $Prefix.screen.png ($($screen.Width) x $($screen.Height))"
}

Get-Content -Path $out
exit 0
