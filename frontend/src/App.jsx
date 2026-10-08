import React, { useState, useEffect } from 'react';
import Header from './components/Header';
import Home from './pages/Home';
import Run from './pages/Run';
import Results from './pages/Results';

const App = () => {
  const [currentPage, setCurrentPage] = useState('home');
  const [currentRunId, setCurrentRunId] = useState('');

  const navigate = (page) => {
    window.scrollTo({ top: 0, left: 0, behavior: 'instant' });
    setCurrentPage(page);
  };

  useEffect(() => {
    window.scrollTo({ top: 0, left: 0, behavior: 'instant' });
  }, [currentPage]);

  return (
    <div className="main-container">
      <Header currentPage={currentPage} onNavigate={navigate} />
      <div className="content-container">
        {currentPage === 'home' && <Home onNavigate={navigate} />}
        {currentPage === 'run' && <Run onNavigate={navigate} setCurrentRunId={setCurrentRunId} />}
        {currentPage === 'results' && <Results currentRunId={currentRunId} />}
      </div>
    </div>
  );
};

export default App;
