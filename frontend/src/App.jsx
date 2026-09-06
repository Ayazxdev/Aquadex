import React, { useState } from 'react';
import Header from './components/Header';
import Home from './pages/Home';
import Run from './pages/Run';
import Results from './pages/Results';

const App = () => {
  const [currentPage, setCurrentPage] = useState('home');
  const [currentRunId, setCurrentRunId] = useState('');

  const navigate = (page) => {
    setCurrentPage(page);
  };

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
