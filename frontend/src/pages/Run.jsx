import React, { useState, useRef } from 'react';
import { Upload, FileText, CheckCircle, AlertCircle, Loader, RefreshCw } from 'lucide-react';
import { API_BASE } from '../utils/config';
import './run.css';

const Run = ({ onNavigate, setCurrentRunId }) => {
  const [files, setFiles] = useState([]);
  const [marker, setMarker] = useState('16S');
  const [readType, setReadType] = useState('short');
  const [isUploading, setIsUploading] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [status, setStatus] = useState('');
  const [runId, setRunId] = useState('');
  const [dragActive, setDragActive] = useState(false);

  // Use a ref for the active runId so pollStatus never captures a stale closure
  const activeRunIdRef = useRef('');
  const pollTimerRef = useRef(null);

  const resetAll = () => {
    if (pollTimerRef.current) clearTimeout(pollTimerRef.current);
    setFiles([]);
    setIsUploading(false);
    setIsRunning(false);
    setStatus('');
    setRunId('');
    activeRunIdRef.current = '';
  };

  const handleFileChange = (e) => {
    setFiles(Array.from(e.target.files));
  };

  const handleDrag = (e) => {
    e.preventDefault();
    e.stopPropagation();
    if (e.type === 'dragenter' || e.type === 'dragover') {
      setDragActive(true);
    } else if (e.type === 'dragleave') {
      setDragActive(false);
    }
  };

  const handleDrop = (e) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      setFiles(Array.from(e.dataTransfer.files));
    }
  };

  const removeFile = (index) => {
    setFiles(files.filter((_, i) => i !== index));
  };

  const uploadFiles = async () => {
    if (files.length === 0) return;

    setIsUploading(true);
    setStatus('Uploading files…');
    const formData = new FormData();
    files.forEach(file => formData.append('files', file));

    try {
      const response = await fetch(`${API_BASE}/upload`, {
        method: 'POST',
        body: formData,
      });
      if (!response.ok) {
        const errText = await response.text();
        setStatus(`Upload failed (${response.status}): ${errText}`);
        return;
      }
      const result = await response.json();
      const serverRunId = result.run_id || result.runId || result.runID || result.id;
      if (!serverRunId) {
        setStatus('Upload succeeded but no run_id returned from server.');
        return;
      }
      activeRunIdRef.current = serverRunId;
      setRunId(serverRunId);
      setCurrentRunId(serverRunId);
      setStatus(`Files uploaded: ${(result.saved_files || []).join(', ')}`);
    } catch (error) {
      setStatus(`Upload error: ${error.message}`);
    } finally {
      setIsUploading(false);
    }
  };

  // Accept the resolved runId as a parameter to avoid stale closure
  const pollStatus = async (id, retries = 0) => {
    if (!id) return;
    try {
      const response = await fetch(`${API_BASE}/status/${id}`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const result = await response.json();
      setStatus(`Status: ${result.status} (${Math.round((result.progress || 0) * 100)}%)`);
      if (result.status === 'completed') {
        setIsRunning(false);
        setStatus('Pipeline completed! View results.');
      } else if (result.status === 'failed') {
        setIsRunning(false);
        setStatus(`Pipeline failed: ${result.message || 'Unknown error'}`);
      } else {
        // Continue polling — pass id explicitly so closure is never stale
        pollTimerRef.current = setTimeout(() => pollStatus(id, retries), 5000);
      }
    } catch (error) {
      if (retries < 3) {
        pollTimerRef.current = setTimeout(() => pollStatus(id, retries + 1), 3000);
      } else {
        setStatus(`Status polling error: ${error.message}`);
        setIsRunning(false);
      }
    }
  };

  const startPipeline = async () => {
    const id = activeRunIdRef.current || runId;
    if (!id) return;

    setIsRunning(true);
    setStatus('Launching pipeline…');
    try {
      const response = await fetch(`${API_BASE}/run`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          run_id: id,
          marker: marker,
          read_type: readType,
          options: {}
        }),
      });
      if (!response.ok) {
        const errText = await response.text();
        setStatus(`Failed to start pipeline (${response.status}): ${errText}`);
        setIsRunning(false);
        return;
      }
      setStatus('Pipeline started! Monitoring progress…');
      pollStatus(id);
    } catch (error) {
      setStatus(`Pipeline error: ${error.message}`);
      setIsRunning(false);
    }
  };

  const getStatusIcon = () => {
    if (isUploading || isRunning) return <Loader className="status-icon spinning" />;
    if (status.includes('completed')) return <CheckCircle className="status-icon success" />;
    if (status.toLowerCase().includes('error') || status.includes('failed')) return <AlertCircle className="status-icon error" />;
    return null;
  };

  const isFailed = status.toLowerCase().includes('failed') || status.toLowerCase().includes('error');

  return (
    <div className="run-page">
      <div className="run-container">
        <div className="run-header">
          <h1 className="run-title">Run eDNA Analysis</h1>
          <p className="run-subtitle">Upload your sequencing files and configure your analysis pipeline</p>
        </div>

        <div className="run-content">
          {/* Upload Section */}
          <div className="upload-section">
            <h2 className="section-title">Upload Files</h2>
            <div
              className={`upload-area ${dragActive ? 'drag-active' : ''}`}
              onDragEnter={handleDrag}
              onDragLeave={handleDrag}
              onDragOver={handleDrag}
              onDrop={handleDrop}
            >
              <input
                type="file"
                multiple
                accept=".fasta,.fa,.fas,.fastq,.fq,.fasta.gz,.fa.gz,.fas.gz,.fastq.gz,.fq.gz"
                className="file-input"
                id="file-upload"
                onChange={handleFileChange}
              />
              <label htmlFor="file-upload" className="upload-label">
                <Upload className="upload-icon" />
                <span className="upload-text">Drag &amp; drop your sequencing files here</span>
                <span className="upload-subtext">or click to browse</span>
                <span className="upload-formats">Supported: .fasta, .fa, .fastq, .fq (including .gz compressed)</span>
              </label>
            </div>

            {files.length > 0 && (
              <div className="files-list">
                <h3 className="files-list-title">Selected Files ({files.length})</h3>
                {files.map((file, index) => (
                  <div key={index} className="file-item">
                    <FileText className="file-icon" />
                    <span className="file-name">{file.name}</span>
                    <span className="file-size">{(file.size / 1024).toFixed(1)} KB</span>
                    <button
                      className="file-remove"
                      onClick={() => removeFile(index)}
                      aria-label="Remove file"
                    >
                      ×
                    </button>
                  </div>
                ))}
              </div>
            )}
          </div>

          {/* Configuration Section */}
          <div className="config-section">
            <h2 className="section-title">Configuration</h2>
            <div className="config-grid">
              <div className="config-group">
                <label className="config-label" htmlFor="marker">Marker Gene</label>
                <select
                  id="marker"
                  className="config-select"
                  value={marker}
                  onChange={(e) => setMarker(e.target.value)}
                >
                  <option value="16S">16S rRNA</option>
                  <option value="18S">18S rRNA</option>
                  <option value="COI">COI</option>
                </select>
              </div>
              <div className="config-group">
                <label className="config-label" htmlFor="readType">Read Type</label>
                <select
                  id="readType"
                  className="config-select"
                  value={readType}
                  onChange={(e) => setReadType(e.target.value)}
                >
                  <option value="short">Short reads (Illumina)</option>
                  <option value="long">Long reads (Nanopore/PacBio)</option>
                </select>
              </div>
            </div>
          </div>

          {/* Status Section */}
          {status && (
            <div className="status-section">
              <div className="status-content">
                {getStatusIcon()}
                <span className="status-text">{status}</span>
              </div>
            </div>
          )}

          {/* Action Buttons */}
          <div className="action-section">
            {!runId ? (
              <button
                className={`btn-primary ${isUploading ? 'loading' : ''}`}
                onClick={uploadFiles}
                disabled={isUploading || files.length === 0}
              >
                {isUploading ? (
                  <>
                    <Loader className="btn-icon spinning" />
                    Uploading...
                  </>
                ) : (
                  'Upload Files'
                )}
              </button>
            ) : status.includes('completed') ? (
              <div style={{ display: 'flex', gap: '12px', flexWrap: 'wrap' }}>
                <button
                  className="btn-primary"
                  onClick={() => {
                    setCurrentRunId(runId);
                    onNavigate('results');
                  }}
                >
                  View Results
                </button>
                <button className="btn-secondary" onClick={resetAll}>
                  <RefreshCw className="btn-icon" style={{ width: 16, height: 16, marginRight: 6 }} />
                  New Analysis
                </button>
              </div>
            ) : isFailed ? (
              <div style={{ display: 'flex', gap: '12px', flexWrap: 'wrap' }}>
                <button className="btn-primary" onClick={startPipeline} disabled={isRunning}>
                  Retry Pipeline
                </button>
                <button className="btn-secondary" onClick={resetAll}>
                  <RefreshCw className="btn-icon" style={{ width: 16, height: 16, marginRight: 6 }} />
                  Start Over
                </button>
              </div>
            ) : (
              <button
                className={`btn-primary ${isRunning ? 'loading' : ''}`}
                onClick={startPipeline}
                disabled={isRunning}
              >
                {isRunning ? (
                  <>
                    <Loader className="btn-icon spinning" />
                    Running Analysis...
                  </>
                ) : (
                  'Start Pipeline'
                )}
              </button>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default Run;